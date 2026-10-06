"""휠체어 결합형 기립 보조 기구 - 준정적 계산 (표준 라이브러리만 사용)
기준: 기립보조 시트 3D 시뮬레이터 v3 의 구성(계측 스프링 - 랙·피니언·드럼 - 케이블 - 수레 - 보조 스프링 - 지렛대 - 들어올림 링크 - 좌판)
좌표: z = 전방(+), y = 위(+), 단위 m (힘 계산은 N, mm)
실행: python3 휠체어_기구_계산.py [앞 피벗-하중점 거리 b(mm), 기본 250]  (코드 안 변수명은 d)
"""
import math, sys

g0 = 9.81
# ---- v3 에서 그대로 가져온 값 ----
SEAT_FRAC = 0.70      # 좌판이 받는 체중 비율
L = 380.0             # 앞 피벗 - 들어올림 링크 [mm]
a = 350.0             # 지렛대 길이 [mm]
F_RET = 3.0           # 복귀 스프링 힘(정하중) [N]
R_PIN = 6.0           # 피니언 반지름 [mm]
H_UP = 132.0          # 좌판을 올렸을 때 링크 연결점 높이 [mm] (기울기 20.3도)
d = float(sys.argv[1]) if len(sys.argv) > 1 else 250.0   # 앞 피벗 - 하중점(좌골결절) [mm]
# ---- 휠체어에 맞춰 다시 잡은 값 ----
F = 450.0             # 보조 스프링 힘 [N] (v3: 200)
i = 4.0               # 증폭비 (드럼 / 피니언) (v3: 8)
RATIO = 0.35          # 설계 보조율 (9/28 회의록: 하지 근력의 약 35 %)
PITCH = 2.5           # 톱니판 피치 [mm] (v3: 5)
G = RATIO / (1 - 2 * RATIO)                 # g = i F / (a k)
k = i * F / (a * G)                         # 계측 스프링 강성 [N/mm]
s0 = a * G * i * F_RET / (F * (1 + G))      # 수레 초기 위치 [mm] (복귀 스프링 힘을 상쇄하는 값)

# ---- 휠체어 안 배치 (WYK874-41, 좌석 레일 높이 0.455 m) ----
P = (0.19, 0.455)     # 앞 피벗 (z, y)
O = (0.17, 0.31)      # 지렛대 피벗
Y_BASE = 0.115        # 보조 스프링 아래 끝 높이 (수레 윗면)
LEVER_T = 0.025       # 지렛대 두께

def S(th): return (P[0] - L / 1000 * math.cos(th), P[1] + L / 1000 * math.sin(th))
def E(ph): return (O[0] - a / 1000 * math.cos(ph), O[1] - a / 1000 * math.sin(ph))
def dist(p, q): return math.hypot(p[0] - q[0], p[1] - q[1])
def bisect(f, lo, hi):
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if (f(lo) > 0) == (f(mid) > 0): lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

TH_UP = math.asin(H_UP / L)
LK = dist(S(TH_UP), E(0.0))                                   # 들어올림 링크 길이
def theta_of(ph): return bisect(lambda th: dist(S(th), E(ph)) - LK, -0.2, 0.6)
PH_MAX = bisect(lambda ph: dist(S(0.0), E(ph)) - LK, 0.0, 0.7)  # 좌판 수평일 때 지렛대 각

def link_load(m): return SEAT_FRAC * m * g0 * d / L           # R [N]
def measure(m):
    R = link_load(m)
    dl = (R - F * s0 / a - i * F_RET) / (k * (1 + 2 * G))     # 계측 스프링 눌림 [mm]
    s = s0 + i * dl                                           # 평형 수레 위치 [mm]
    s_lock = s0 + math.floor((s - s0) / PITCH) * PITCH        # 걸쇠에 잠긴 위치
    return R, dl, s, s_lock
def assist(s, ph):
    """수레 위치 s, 지렛대 각 ph 에서 링크 위치 환산 보조력 [N] (가상일: F * d(접촉점 높이)/d(좌판 각))"""
    h = 1e-5
    dph_dth = -2 * h / (theta_of(ph + h) - theta_of(ph - h))   # 좌판이 올라가면(th 증가) 지렛대 각 ph 는 줄어든다
    M = F * (s / 1000) / math.cos(ph) ** 2 * dph_dth          # 앞 피벗 기준 모멘트 [N·m]
    return M / (L / 1000 * math.cos(theta_of(ph)))

if __name__ == '__main__':
    print(f"g = {G:.4f}, k = {k:.3f} N/mm, s0 = {s0:.2f} mm, 드럼 반지름 = {i * R_PIN:.0f} mm")
    print(f"좌판 올림각 = {math.degrees(TH_UP):.2f} deg, 링크 길이 = {LK * 1000:.1f} mm, 좌판 수평 시 지렛대 각 = {math.degrees(PH_MAX):.2f} deg")
    L_ext = O[1] - LEVER_T / 2 - Y_BASE
    print(f"보조 스프링 늘인 길이(지렛대 수평) = {L_ext * 1000:.1f} mm")
    print("체중 | R[N] | 눌림[mm] | s평형 | s잠김 | 보조력(계측) | 보조율(계측) | 보조력(좌판 수평) | 보조율(수평) | 스프링 행정[mm] | 걸쇠 하중[N]")
    for m in range(40, 121, 10):
        R, dl, s, sl = measure(m)
        ph_m = math.asin(dl / a)
        A1 = assist(sl, ph_m); A2 = assist(sl, PH_MAX)
        stroke = sl * math.tan(PH_MAX)
        print(f"{m:4d} | {R:6.1f} | {dl:5.1f} | {s:6.1f} | {sl:6.1f} | {A1:6.1f} | {A1 / R * 100:5.1f} % | {A2:6.1f} | {A2 / R * 100:5.1f} % | {stroke:5.1f} | {F * math.tan(PH_MAX) + F_RET:5.1f}")
    R, dl, s, sl = measure(120)
    print(f"조건: 2*행정+50 = {2 * sl * math.tan(PH_MAX) + 50:.1f} mm <= 늘인 길이 {L_ext * 1000:.1f} mm ?")
    # 선형식 검산: 보조율 = g / (1 + 2g)
    print(f"선형식 보조율 g/(1+2g) = {G / (1 + 2 * G) * 100:.2f} %, F*s/a /R (80 kg) = {F * measure(80)[2] / a / link_load(80) * 100:.2f} %")
