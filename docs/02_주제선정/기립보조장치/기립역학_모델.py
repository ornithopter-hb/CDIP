"""기립(sit-to-stand) 2D 시상면 4분절 모델 - 최소 관절모멘트 궤적 최적화 + 역동역학
좌표: x = 전방(+), y = 위(+) (3D 휠체어 좌표의 z = x)
분절: 발(지면 고정), 하퇴, 대퇴, HAT. 양쪽 다리 합산(좌우 대칭).
"""
import json, math, sys
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize

g = 9.81
H, M = 1.70, 70.0                               # [가정] 신장 1.70 m, 체중 70 kg
h_a = 0.039 * H                                 # Winter Fig 4.1
L_s = (0.285 - 0.039) * H
L_t = (0.530 - 0.285) * H
L_T = (0.818 - 0.530) * H
L_heel = 0.0454 / 1.74 * H                      # Yoshioka et al. 2007 (H = 1.74 m) 비율
L_toe = 0.2018 / 1.74 * H
m_f, m_s, m_th, m_T = 2 * 0.0145 * M, 2 * 0.0465 * M, 2 * 0.100 * M, 0.678 * M   # Winter Table 4.1
c_s, c_th, c_T = 0.433, 0.433, 0.626
I_s = m_s * (0.302 * L_s) ** 2
I_th = m_th * (0.323 * L_t) ** 2
I_T = m_T * (0.496 * L_T) ** 2

seat_y = 0.47                                   # WYK874-41 좌석 높이
hip0 = np.array([0.06, seat_y + 0.09])          # [가정] 대전자: 좌면 위 9 cm, 좌석 앞끝(z=0.20)에서 14 cm 뒤 (앞으로 당겨 앉음)
a_s0 = math.radians(float(sys.argv[1]) if len(sys.argv) > 1 else 15.0)   # [가정] 초기 하퇴 전방경사(배측굴곡)

cb = (hip0[1] - h_a - L_s * math.cos(a_s0)) / L_t
b0 = -math.acos(cb)
A = np.array([hip0[0] - L_s * math.sin(a_s0) - L_t * math.sin(b0), h_a])
kap0 = a_s0 - b0
heel_x, toe_x = A[0] - L_heel, A[0] + L_toe

T_mv = 1.6                                       # Yoshioka 2014: 좌석 40·50 cm에서 1.6 s
t_so, t_p2, t_p3 = 0.28 * T_mv, 0.46 * T_mv, T_mv   # Schenkman 1990: 28 / 18 / 54 %
t_mid = 0.5 * (t_p2 + t_p3)
t_end = T_mv + 0.5
dt = 0.002
ts = np.arange(0, t_end + 1e-9, dt)
after = ts >= t_so - 1e-9

def kin(x):
    th1, th2, th3, k2, k3, a2, a3 = x
    tr = CubicSpline([0, t_so, t_p2, t_mid, t_p3], np.radians([0, th1, th2, th3, 0]), bc_type='clamped')
    kn = CubicSpline([t_so, t_p2, t_mid, t_p3], [kap0, math.radians(k2), math.radians(k3), 0], bc_type='clamped')
    an = CubicSpline([t_so, t_p2, t_mid, t_p3], [a_s0, math.radians(a2), math.radians(a3), 0], bc_type='clamped')
    tc = np.clip(ts, 0, t_p3)
    gm = tr(tc); gd = np.where(ts <= t_p3, tr(tc, 1), 0); gdd = np.where(ts <= t_p3, tr(tc, 2), 0)
    tl = np.clip(ts, t_so, t_p3); on = (ts >= t_so - 1e-9) & (ts <= t_p3)
    kap = kn(tl); kd = np.where(on, kn(tl, 1), 0); kdd = np.where(on, kn(tl, 2), 0)
    a = an(tl); ad = np.where(on, an(tl, 1), 0); add = np.where(on, an(tl, 2), 0)
    return gm, gd, gdd, kap, kd, kdd, a, ad, add

def unit(th): return np.stack([np.sin(th), np.cos(th)], -1)
def dunit(th): return np.stack([np.cos(th), -np.sin(th)], -1)

def dynamics(x):
    gm, gd, gdd, kap, kd, kdd, a, ad, add = kin(x)
    b, bd, bdd = a - kap, ad - kd, add - kdd
    # 위치
    K = A + L_s * unit(a); P = K + L_t * unit(b); S = P + L_T * unit(gm)
    Cs = A + (1 - c_s) * L_s * unit(a); Cth = K + (1 - c_th) * L_t * unit(b); CT = P + c_T * L_T * unit(gm)
    # 가속도 (해석적): r = L u(θ) → r'' = L(θ'' du - θ'^2 u)
    def acc(L, th, thd, thdd): return L * (thdd[:, None] * dunit(th) - (thd ** 2)[:, None] * unit(th))
    aK = acc(L_s, a, ad, add); aP = aK + acc(L_t, b, bd, bdd)
    aCs = acc((1 - c_s) * L_s, a, ad, add); aCth = aK + acc((1 - c_th) * L_t, b, bd, bdd)
    aCT = aP + acc(c_T * L_T, gm, gd, gdd)
    cross = lambda r, F: r[:, 0] * F[:, 1] - r[:, 1] * F[:, 0]
    gv = np.array([0.0, g])
    wdT, wdth, wds = -gdd, -bdd, -add                     # 반시계 + 각가속도
    F_h = m_T * (aCT + gv)
    M_h = I_T * wdT - cross(P - CT, F_h)                  # 고관절 신전 +
    F_k = m_th * (aCth + gv) + F_h
    M_k = I_th * wdth + M_h + cross(P - Cth, F_h) - cross(K - Cth, F_k)
    F_a = m_s * (aCs + gv) + F_k
    M_a = I_s * wds + M_k + cross(K - Cs, F_k) - cross(A[None, :] - Cs, F_a)
    G = F_a + np.array([0.0, m_f * g])
    Cf = np.array([heel_x + 0.5 * (L_heel + L_toe), 0.5 * h_a])
    xcop = A[0] + (M_a + (Cf[0] - A[0]) * m_f * g - A[1] * G[:, 0]) / G[:, 1]
    com = (m_f * Cf + m_s * Cs + m_th * Cth + m_T * CT) / M
    vcom = np.gradient(com, dt, axis=0)
    return dict(gm=gm, kap=kap, a=a, K=K, P=P, S=S, CT=CT, com=com, vcom=vcom, G=G, xcop=xcop,
                hip=M_h / M, knee=-M_k / M, ankle=M_a / M, Fz=G[:, 1] / (M * g))

def cost(x):
    r = dynamics(x)
    J = np.sum(r['hip'] ** 2) + np.sum(r['knee'][after] ** 2) + np.sum(r['ankle'][after] ** 2)
    return J * dt

margin = 0.01
def cons(x):
    r = dynamics(x)
    c = r['xcop'][after]
    return np.concatenate([c - (heel_x + margin), (toe_x - margin) - c, r['G'][after, 1] - 50.0,
                           np.radians(np.array([x[3], x[4]])) - 0.0,
                           np.radians(TMAX) - r['gm']])

x0 = np.array([22, 30, 15, math.degrees(kap0) - 15, 40, math.degrees(a_s0) + 5, 12])
TMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 35.0   # 몸통 최대 경사 상한 (Yoshioka 2014 참고)
bnds = [(0, TMAX), (0, TMAX), (-5, TMAX), (0, 130), (0, 130), (math.degrees(a_s0), 45), (-10, 45)]
sol = minimize(cost, x0, method='SLSQP', bounds=bnds,
               constraints=[{'type': 'ineq', 'fun': cons}], options=dict(maxiter=400, ftol=1e-8))
x = sol.x
r = dynamics(x)
i_so = int(round(t_so / dt)); i2 = int(round(t_p2 / dt))
w0 = math.sqrt(g / 0.9)
def pk(key, sl=slice(i_so, None)):
    v = r[key][sl]; i = np.argmax(v); return round(float(v[i]), 3), round(float(ts[sl][i]), 3)
out = dict(
    a_s0=math.degrees(a_s0), kap0=round(math.degrees(kap0), 1), success=bool(sol.success), msg=sol.message,
    x=np.round(x, 2).tolist(), cost=round(float(sol.fun), 3),
    ankle_x=float(A[0]), heel_x=float(heel_x), toe_x=float(toe_x), hip0=hip0.tolist(),
    t_so=t_so, t_p2=t_p2, t_p3=t_p3,
    trunk_max=round(float(np.degrees(r['gm'].max())), 1), t_trunk_max=round(float(ts[np.argmax(r['gm'])]), 3),
    ankle_df_max=round(float(np.degrees(r['a'].max())), 1), t_df_max=round(float(ts[np.argmax(r['a'])]), 3),
    Fz=pk('Fz'), hip=pk('hip', slice(None)), knee=pk('knee'), ankle=pk('ankle'),
    ankle_min=round(float(r['ankle'][after].min()), 3),
    cop_range=[round(float(r['xcop'][after].min() - A[0]), 3), round(float(r['xcop'][after].max() - A[0]), 3)],
    seatoff=dict(trunk=round(float(np.degrees(r['gm'][i_so])), 1),
                 trunk_rate=round(float(np.degrees(np.gradient(r['gm'], dt)[i_so])), 1),
                 com_rel_ankle=round(float(r['com'][i_so, 0] - A[0]), 3),
                 com_rel_heel=round(float(r['com'][i_so, 0] - heel_x), 3),
                 vx=round(float(r['vcom'][i_so, 0]), 3), vy=round(float(r['vcom'][i_so, 1]), 3),
                 xcom_rel_ankle=round(float(r['com'][i_so, 0] + r['vcom'][i_so, 0] / w0 - A[0]), 3),
                 hip=round(float(r['hip'][i_so]), 3), knee=round(float(r['knee'][i_so]), 3),
                 ankle=round(float(r['ankle'][i_so]), 3), Fz=round(float(r['Fz'][i_so]), 3),
                 cop_rel_ankle=round(float(r['xcop'][i_so] - A[0]), 3)),
    p2end=dict(com_rel_ankle=round(float(r['com'][i2, 0] - A[0]), 3)),
    com_start=np.round(r['com'][0], 3).tolist(), com_end=np.round(r['com'][-1], 3).tolist(),
    vx_max=round(float(r['vcom'][:, 0].max()), 3), vy_max=round(float(r['vcom'][:, 1].max()), 3),
    t_vy_max=round(float(ts[np.argmax(r['vcom'][:, 1])]), 3),
)
print(json.dumps(out, ensure_ascii=False))
if len(sys.argv) > 2 and sys.argv[2] != '-':
    step = 5
    def L(v): return [None if not np.isfinite(z) else round(float(z), 4) for z in v]
    ser = dict(t=np.round(ts[::step], 3).tolist(), A=np.round(A, 4).tolist(), heel=float(heel_x), toe=float(toe_x),
               K=np.round(r['K'][::step], 4).tolist(), P=np.round(r['P'][::step], 4).tolist(),
               S=np.round(r['S'][::step], 4).tolist(), com=np.round(r['com'][::step], 4).tolist(),
               Fz=L(np.where(after, r['Fz'], np.nan)[::step]),
               Fx=L(np.where(after, r['G'][:, 0] / (M * g), np.nan)[::step]),
               cop=L(np.where(after, r['xcop'], np.nan)[::step]),
               hip=L(r['hip'][::step]), knee=L(np.where(after, r['knee'], np.nan)[::step]),
               ankle=L(np.where(after, r['ankle'], np.nan)[::step]),
               phases=[t_so, t_p2, t_p3], summary=out)
    json.dump(ser, open(sys.argv[2], 'w'), ensure_ascii=False)

# ---- 추가 검산: 전체 운동량 변화율 = 외력 (이탈 후, 전 해상도) ----
acom = np.gradient(np.gradient(r['com'], dt, axis=0), dt, axis=0)
Gt = M * (acom + np.array([0.0, g]))
mask = after.copy()
for tb in (t_so, t_p2, t_mid, t_p3):
    mask &= np.abs(ts - tb) > 0.01
mask[:3] = False; mask[-3:] = False
chk = float(np.max(np.abs(Gt[mask] - r['G'][mask])) / (M * g))
# 준정적 이탈: 다리 초기자세 고정, 몸통각만 변할 때 CoM이 뒤꿈치·발목 위에 오는 몸통각
def com_static(gdeg):
    gmm = math.radians(gdeg)
    K = A + L_s * np.array([math.sin(a_s0), math.cos(a_s0)]); P = hip0
    S_ = P + L_T * np.array([math.sin(gmm), math.cos(gmm)])
    CT_ = P + c_T * (S_ - P); Cth_ = P + c_th * (K - P); Cs_ = K + (1 - c_s) * (A - K)
    Cf_ = np.array([heel_x + 0.5 * (L_heel + L_toe), 0.5 * h_a])
    return ((m_f * Cf_ + m_s * Cs_ + m_th * Cth_ + m_T * CT_) / M)[0]
def solve_trunk(target):
    lo, hi = 0.0, 90.0
    if com_static(hi) < target: return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if com_static(mid) > target: hi = mid
        else: lo = mid
    return round(0.5 * (lo + hi), 1)
extra = dict(newton_check_BW=round(chk, 5), static_trunk_heel=solve_trunk(heel_x), static_trunk_ankle=solve_trunk(A[0]))
print(json.dumps(extra))
if len(sys.argv) > 2 and sys.argv[2] != '-':
    ser = json.load(open(sys.argv[2])); ser['summary'].update(extra); json.dump(ser, open(sys.argv[2], 'w'), ensure_ascii=False)
