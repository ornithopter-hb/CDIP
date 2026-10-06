"""휠체어 측면 프레임 내장형 기립 보조 기구 - 좌판 최대 33도, 링크식(줄 없음), 레버 없음, 바닥 발(자동 브레이크). 준정적 계산 (표준 라이브러리만 사용)
구조
  지렛대 축 = 뒷바퀴 축 O. 지렛대와 한 몸인 부채꼴 팔의 호 홈을 슬라이더가 탄다. 보조 스프링(가스, 눕힘)은 고정 핀 B 와 슬라이더 사이.
  계측 팔: B 에 단 지렛대. 긴 팔(ELL) 끝이 슬라이더를 밀고, 짧은 팔(RS) 끝에 받침이 선다. 짧은 팔은 계측 범위 중간에서 지렛대 코와 나란하다.
  받침: 짧은 팔 끝에 핀으로 선 막대(LS). 윗끝 롤러가 지렛대 코 아랫면에 닿는다. 코 아랫면의 법선에서 뒤로 PSI0 기울어 서 있다.
        받침은 아래 핀 밑으로 꼬리(HT)가 나와 있고, 지연 댐퍼(양끝 핀: 프레임 DF - 꼬리 핀)가 꼬리를 붙잡는다. 댐퍼는 눌리는 쪽만 저항한다.
        댐퍼는 다 늘어난 길이(L0)에서 멈추므로, 하중이 없을 때 받침은 댐퍼 안의 약한 스프링으로 PSI0 에 서 있다.
        댐퍼가 XL 만큼 눌리면 걸쇠가 물리고, XR 만큼 눌리면 댐퍼가 놓는다(그 뒤 행정은 자유).
        그 뒤 받침은 코 밑에 붙은 채 뒤로 누웠다가, 일어설 때 코를 따라 다시 선다.
  바닥 발: 뒤 가로 파이프에 자유롭게 도는 다리. 좌판이 올라가면 내려와 뒷바퀴 뒤 바닥을 짚는다(뒤로 밀림·뒤로 들림 방지).
           지렛대 축의 크랭크와 다리의 귀를 장공 막대로 이어, 좌판이 내려가 있으면 막대가 다리를 들고 있다.
좌표: z = 전방(+), y = 위(+), 단위 m, N. 힘·강성은 좌우 합.
실행: python3 휠체어_측면내장_33도_계산.py [앞 피벗-하중점 거리 b(mm), 기본 250]
"""
import math, sys
g0=9.81
def dist(p,q): return math.hypot(p[0]-q[0],p[1]-q[1])
def bis(f,lo,hi,n=48):
    flo=f(lo)
    for _ in range(n):
        m=0.5*(lo+hi); fm=f(m)
        if (flo>0)==(fm>0): lo,flo=m,fm
        else: hi=m
    return 0.5*(lo+hi)
class M:
    def __init__(s, TH_UP=math.radians(33), AR=0.24, LP=0.20, EX=None, F=800.0, B=(0.07,0.30), O=(-0.17,0.295), P=(0.19,0.455), DROP=0.012,
                 I=8.0, U0=0.025, PSI0=math.radians(6), PSIL=math.radians(10), PSIR=math.radians(14), HT=0.04, LD=0.11, RR=0.005, HW=0.0135, target=0.35, b=0.25, DMID=None):
        s.__dict__.update(TH_UP=TH_UP,AR=AR,LP=LP,F=F,B=B,O=O,P=P,DROP=DROP,I=I,U0=U0,PSI0=PSI0,PSIL=PSIL,PSIR=PSIR,HT=HT,LD=LD,RR=RR,HW=HW,target=target,b=b)
        s.FRAC=0.7; s.FRET=3.0
        s.LK=dist(s.S(0),s.E(0)); s.PH_UP=bis(lambda ph: s.th_of(ph)-TH_UP,0,1.3)
        s.ELL=dist(B,O); s.BASE=math.atan2(O[1]-B[1],O[0]-B[0]); s.RS=s.ELL/I
        s.DMID=DMID if DMID is not None else 0.25          # 짧은 팔이 코 방향과 나란해지는 회전각 (계측 범위 중간쯤)
        s.A0=s.PH_UP+s.DMID                                # 짧은 팔의 처음 각 (수평에서 위로, 앞쪽)
        # 받침 길이: 계측 팔이 처음 자리(D=0)이고 받침이 PSI0 일 때 지렛대가 PH_UP 가 되도록
        bp=s.Bp(0); dirv=(math.cos(s.PH_UP),math.sin(s.PH_UP))
        d=(bp[0]-O[0])*dirv[1]-(bp[1]-O[1])*dirv[0]       # 짧은 팔 끝에서 지렛대 중심선까지 수직 거리
        s.LS=(d-HW-RR)/math.cos(PSI0)
        # 지연 댐퍼의 프레임 쪽 핀 DF: 받침 각이 PSI0 그대로일 때, 계측 팔이 돌아도(체중이 달라도) 꼬리 핀까지 거리가 가장 덜 변하는 방향을 고른다
        t0=s.tail(0,PSI0); best=None
        for a10 in range(0,500):
            al=math.radians(a10/10); df=(t0[0]+LD*math.cos(al), t0[1]+LD*math.sin(al))
            dl=[dist(s.tail(D,PSI0),df)-LD for D in [i*0.03 for i in range(1,16)]]
            err=max(abs(v) for v in dl)
            if best is None or err<best[0]: best=(err,al,df,min(dl),max(dl))
        s.DFERR,s.DFANG,s.DF,s.DFMIN,s.DFMAX=best; s.L0=LD
        s.XL=s.L0-dist(s.tail(0,PSIL),s.DF); s.XR=s.L0-dist(s.tail(0,PSIR),s.DF)
    def S(s,th): return (s.P[0]-s.LP*math.cos(th), s.P[1]+s.LP*math.sin(th)-s.DROP)
    def E(s,ph): return (s.O[0]+s.AR*math.cos(ph), s.O[1]+s.AR*math.sin(ph))
    def th_of(s,ph): return bis(lambda th: dist(s.S(th),s.E(ph))-s.LK,-0.15,0.9,42)
    def Q0(s,u):
        be=s.BASE-u/s.ELL; return (s.B[0]+s.ELL*math.cos(be), s.B[1]+s.ELL*math.sin(be))
    def Q(s,u,ph):
        q=s.Q0(u); d=ph-s.PH_UP; x=q[0]-s.O[0]; y=q[1]-s.O[1]
        return (s.O[0]+x*math.cos(d)-y*math.sin(d), s.O[1]+x*math.sin(d)+y*math.cos(d))
    def Ls(s,u,ph): return dist(s.Q(u,ph),s.B)
    def Ms(s,u,ph):
        q=s.Q(u,ph); L=dist(q,s.B); fx=s.F*(q[0]-s.B[0])/L; fy=s.F*(q[1]-s.B[1])/L
        return (q[0]-s.O[0])*fy-(q[1]-s.O[1])*fx
    def Ft(s,u,ph): h=1e-6; return s.F*(s.Ls(u+h,ph)-s.Ls(u-h,ph))/(2*h)
    def Mload(s,m,ph,b=None):
        b=b or s.b; h=1e-4; dth=(s.th_of(ph+h)-s.th_of(ph-h))/(2*h)
        return s.FRAC*m*g0*b*math.cos(s.th_of(ph))*dth
    def Bp(s,D):
        a=s.A0-D; return (s.B[0]+s.RS*math.cos(a), s.B[1]+s.RS*math.sin(a))
    def lever_of(s,D,psi):
        """계측 팔 회전 D, 받침 각 psi(코 법선에서 뒤로) 일 때 지렛대 각과 롤러 중심"""
        bp=s.Bp(D)
        def f(ph):
            c=(bp[0]-s.LS*math.sin(ph+psi), bp[1]+s.LS*math.cos(ph+psi))
            return (c[0]-s.O[0])*math.sin(ph)-(c[1]-s.O[1])*math.cos(ph)-(s.HW+s.RR)
        ph=bis(f,-0.3,s.PH_UP+0.3)
        c=(bp[0]-s.LS*math.sin(ph+psi), bp[1]+s.LS*math.cos(ph+psi))
        return ph,c
    def tail(s,D,psi):
        """받침 꼬리 핀 (아래 핀에서 받침 반대쪽으로 HT)"""
        ph,_=s.lever_of(D,psi); bp=s.Bp(D)
        return (bp[0]+s.HT*math.sin(ph+psi), bp[1]-s.HT*math.cos(ph+psi))
    def psi_at(s,D,x):
        """계측 팔이 D 만큼 돌았을 때 댐퍼가 x 만큼 눌리는 받침 각"""
        return bis(lambda psi: dist(s.tail(D,psi),s.DF)-(s.L0-x), 0.02, 1.2)
    def psi_down(s,ph,D=0.0):
        """댐퍼가 놓은 뒤: 롤러가 코 밑에 붙어 있을 때의 받침 각"""
        bp=s.Bp(D); d=(bp[0]-s.O[0])*math.sin(ph)-(bp[1]-s.O[1])*math.cos(ph)
        return math.acos(max(-1,min(1,(d-s.HW-s.RR)/s.LS)))
    def measure(s,m,kth,tau0,b=None,x=None):
        """걸쇠가 물리는 순간(댐퍼가 XL 만큼 눌림)의 평형. 댐퍼 반력이 계측 팔에 주는 토크를 포함한다."""
        x=s.XL if x is None else x
        def parts(D):
            psi=s.psi_at(D,x); ph,c=s.lever_of(D,psi); u=s.U0+s.ELL*D
            rc=(c[0]-s.O[0])*math.cos(ph)+(c[1]-s.O[1])*math.sin(ph)
            N=(s.Mload(m,ph,b)-s.Ms(u,ph))/rc
            tp=s.tail(D,psi); L=dist(tp,s.DF); uz=(tp[0]-s.DF[0])/L; uy=(tp[1]-s.DF[1])/L      # 댐퍼가 꼬리 핀을 미는 방향
            br=-(math.sin(ph+psi)*uy+math.cos(ph+psi)*uz)                                      # 댐퍼 축이 받침에 수직인 정도 (1 이면 수직)
            Fd=N*s.LS*math.sin(psi)/(s.HT*br)
            Fz=N*math.sin(ph)+Fd*uz; Fy=-N*math.cos(ph)+Fd*uy                                  # 받침 아래 핀이 짧은 팔 끝에 주는 힘
            a=s.A0-D; tl=-(s.RS*math.cos(a)*Fy-s.RS*math.sin(a)*Fz)                            # 계측 팔을 돌리는 토크 (+)
            T=s.FRET-s.Ft(u,ph)
            return tl-(tau0+kth*D+s.ELL*T), ph,u,N,Fd,T,rc,psi,br
        if parts(0.0)[0]<=0: D=0.0
        else: D=bis(lambda D: parts(D)[0],0.0,0.9,44)
        r=parts(D); return dict(D=D,ph=r[1],u=r[2],N=r[3],Fd=r[4],T=r[5],rc=r[6],psi=r[7],br=r[8])
    def er(s,u,m,b=None):
        b=b or s.b
        return s.F*(s.Ls(u,s.PH_UP)-s.Ls(u,0.0))/(s.FRAC*m*g0*b*math.sin(s.TH_UP))
    def fit(s):
        def t0_for(k): return bis(lambda t0: s.er(s.measure(40,k,t0)['u'],40)-s.target,-5.0,30.0,34)
        k=bis(lambda k: s.er(s.measure(120,k,t0_for(k))['u'],120)-s.target,0.5,400.0,34)
        return k,t0_for(k)
    def report(s,name):
        k,t0=s.fit()
        print(f'[{name}] LK={s.LK*1000:.0f} PH_UP={math.degrees(s.PH_UP):.1f} ELL={s.ELL*1000:.0f} RS={s.RS*1000:.1f} LS={s.LS*1000:.1f} kth={k:.2f} Nm/rad tau0={t0:.2f} Nm')
        for m in (40,60,80,100,120):
            r=s.measure(m,k,t0); u=r['u']; q=s.Q(u,s.PH_UP)
            print(f"  {m}: D={math.degrees(r['D']):.1f}deg u={u*1000:.1f} ph={math.degrees(r['ph']):.1f} th={math.degrees(s.th_of(r['ph'])):.1f} N={r['N']:.0f} Fd={r['Fd']:.0f} T={r['T']:.1f} rc={r['rc']*1000:.0f} avg={s.er(u,m)*100:.1f} start={s.Ms(u,0)/s.Mload(m,0)*100:.1f} end={s.Ms(u,s.PH_UP)/s.Mload(m,s.PH_UP)*100:.1f} stroke={(s.Ls(u,s.PH_UP)-s.Ls(u,0))*1000:.1f} Qy={q[1]:.3f} latch={-s.Ft(u,0):.0f} link={s.Mload(m,0)/s.AR:.0f}")
        return k,t0

# ---- 바닥 발 (스프래그 다리) ----
PF = (-0.205, 0.26)                 # 다리 피벗 = 뒤 가로 파이프
LAM = math.radians(25.0)            # 바닥을 짚었을 때 다리가 수직에서 뒤로 기운 각
LF = PF[1] / math.cos(LAM)          # 다리 길이
ZFOOT = PF[0] - PF[1] * math.tan(LAM)
ZR, ZF, WC, ZCG = -0.17, 0.23, 20 * g0, 0.0     # 뒷바퀴·앞바퀴 접지점, 휠체어 무게, 무게중심(가정)
def foot_state(th, m, mus, b=0.25, P=(0.19, 0.455)):
    """좌판 각 th 에서 체중을 다 실었을 때. 다리는 양끝 핀이라 힘이 다리 축 방향으로만 전달된다."""
    V = 0.7 * m * g0; H = max(0.0, V * math.tan(th - math.atan(mus)))
    zs = P[0] - b * math.cos(th); hs = 0.47 + b * math.sin(th)
    Fv = H / math.tan(LAM)                                     # 다리가 받는 수직력 (수평력 H 를 버티려면)
    # 앞바퀴 접지점 기준 모멘트: 뒷바퀴가 뜨지 않는 여유 (+ 이면 뒷바퀴가 바닥에 남음)
    rear = (ZF - zs) * V + (ZF - ZCG) * WC + hs * H - (ZF - PF[0]) * Fv - PF[1] * H
    # 발 접지점 기준 모멘트: 앞바퀴가 뜨지 않는 여유 (뒤로 들림)
    front = (zs - ZFOOT) * V + (ZCG - ZFOOT) * WC - hs * H
    return dict(H=H, Fv=Fv, rear=rear, front=front)
# 바닥 발을 들고 내리는 링크: 지렛대 축의 크랭크(RC) - 장공 막대 - 다리의 귀(RL)
OAX = (-0.17, 0.295); RC, GC0, RL, LAM_STOW = 0.052, math.radians(-185.0), 0.018, math.radians(55.0)
def crank_pin(ph): return (OAX[0] + RC * math.cos(GC0 + ph), OAX[1] + RC * math.sin(GC0 + ph))
def lug(lam): return (PF[0] - RL * math.sin(lam), PF[1] - RL * math.cos(lam))
LROD = dist(crank_pin(0.0), lug(LAM_STOW))
def leg_angle(ph):
    """지렛대 각 ph 에서 다리가 수직에서 뒤로 기운 각. 막대는 당기는 쪽으로만 힘을 전한다(장공)."""
    if dist(crank_pin(ph), lug(LAM)) <= LROD: return LAM            # 발이 바닥에 닿음, 핀은 장공 안에서 미끄러진다
    return bis(lambda lam: dist(crank_pin(ph), lug(lam)) - LROD, LAM, LAM_STOW + 0.2)
def limit(f):
    for d10 in range(0, 900):
        if not f(math.radians(d10 / 10)): return (d10 - 1) / 10
    return 90.0

if __name__ == '__main__':
    b = (float(sys.argv[1]) if len(sys.argv) > 1 else 250.0) / 1000
    s = M(); kth, tau0 = s.fit()          # 설계값은 b = 250 mm 기준
    print(f"좌판 최대각 {math.degrees(s.TH_UP):.0f} deg, 지렛대 회전 {math.degrees(s.PH_UP):.1f} deg, 링크 {s.LK*1000:.0f} mm, 호 반지름(긴 팔) {s.ELL*1000:.0f} mm, 짧은 팔 {s.RS*1000:.0f} mm, 받침 {s.LS*1000:.1f} mm")
    print(f"계측 스프링(비틀림) {kth:.2f} N m/rad, 예압 {tau0:.2f} N m  | raw {kth!r} {tau0!r} {s.PH_UP!r} {s.LS!r} {s.A0!r}")
    r120 = s.measure(120, kth, tau0)
    lt = lambda a: math.log(math.tan(a / 2))
    c_r = 1.0 * r120['N'] * s.LS / (lt(r120['psi']) - lt(s.psi_at(r120['D'], 0.0)))     # 120 kg 에서 걸쇠까지 1.0 초
    print(f"지연 댐퍼: 프레임 쪽 핀 DF = ({s.DF[0]:.3f}, {s.DF[1]:.3f}), 수평에서 {math.degrees(s.DFANG):.1f} deg 위, 길이 {s.L0*1000:.0f} mm, 꼬리 {s.HT*1000:.0f} mm, 받침 각이 그대로일 때 계측 팔 회전에 따른 거리 변화 {s.DFMIN*1000:+.2f} ~ {s.DFMAX*1000:+.2f} mm")
    print(f"  걸쇠 XL = {s.XL*1000:.2f} mm 눌림, 놓음 XR = {s.XR*1000:.2f} mm 눌림, c_r = {c_r:.1f} N m s/rad (댐퍼 축 기준 {c_r/(s.HT*r120['br'])**2/1000:.1f} N s/mm) | raw {c_r!r} DF {s.DF!r}")
    h = 1e-4; dth0 = (s.th_of(h) - s.th_of(-h)) / (2 * h)
    print("체중 | 계측 팔 회전 | 슬라이더 평형/잠김 [mm] | 좌판 각(계측 중) | 평균 | 시작 | 끝 | 스프링 행정 | 받침 N | 댐퍼 힘 | 미는 힘 | 받침 각 처음/걸쇠/놓음 | 걸쇠까지 [s] | 놓을 때까지 [s] | 걸쇠 | 좌판 뒤끝 걸쇠 | 링크")
    for m in range(40, 121, 10):
        r = s.measure(m, kth, tau0, b); u = r['u']; ul = s.U0 + math.floor((u - s.U0) / 0.0025 + 1e-9) * 0.0025
        p0 = s.psi_at(r['D'], 0.0); pL = r['psi']; pR = s.psi_at(r['D'], s.XR)
        tL = c_r / (r['N'] * s.LS) * (lt(pL) - lt(p0)); tR = c_r / (r['N'] * s.LS) * (lt(pR) - lt(p0))
        print(f"{m:4d} | {math.degrees(r['D']):5.1f} | {u*1000:6.1f} / {ul*1000:6.1f} | {math.degrees(s.th_of(r['ph'])):5.1f} | {s.er(ul,m,b)*100:5.1f} % | {s.Ms(ul,0)/s.Mload(m,0,b)*100:5.1f} % | {s.Ms(ul,s.PH_UP)/s.Mload(m,s.PH_UP,b)*100:5.1f} % |"
              f" {(s.Ls(ul,s.PH_UP)-s.Ls(ul,0))*1000:5.1f} | {r['N']:5.0f} | {r['Fd']:5.0f} | {r['T']:5.1f} | {math.degrees(p0):4.1f}/{math.degrees(pL):4.1f}/{math.degrees(pR):4.1f} | {tL:5.2f} | {tR:5.2f} | {-s.Ft(ul,0):5.0f} | {s.Ms(ul,0)/dth0/0.39:5.0f} | {s.Mload(m,0,b)/s.AR:5.0f}")
    pd = s.psi_down(0.0); bp = s.Bp(0); cd = (bp[0]-s.LS*math.sin(pd), bp[1]+s.LS*math.cos(pd)); td = (bp[0]+s.HT*math.sin(pd), bp[1]-s.HT*math.cos(pd))
    print(f"좌판 수평: 받침 각 {math.degrees(pd):.0f} deg, 롤러 중심 ({cd[0]:.3f}, {cd[1]:.3f}) 아랫면 {cd[1]-s.RR:.4f} m (아래 레일 윗면 0.271), 꼬리 핀 ({td[0]:.3f}, {td[1]:.3f}), 댐퍼 길이 {dist(td,s.DF)*1000:.0f} mm")
    lens=[dist((s.Bp(0)[0]+s.HT*math.sin(ph+s.psi_down(ph)), s.Bp(0)[1]-s.HT*math.cos(ph+s.psi_down(ph))),s.DF) for ph in [s.PH_UP*i/60 for i in range(61)]]
    print(f"  댐퍼 길이 범위 {min(lens)*1000:.0f}-{max(lens)*1000:.0f} mm")
    # 좌판 뒤끝 윗모서리와 등받이 바 사이
    gap=min(math.hypot(0.19-0.39*math.cos(th)-0.015*math.sin(th)-(-0.194), 0.455+0.39*math.sin(th)+0.015*math.cos(th)-0.625) for th in [s.TH_UP*i/200 for i in range(201)])-0.009
    print(f"좌판 뒤끝 윗모서리 - 등받이 바 최소 간격 {gap*1000:.0f} mm")
    q = s.Q(s.measure(120, kth, tau0)['u'], s.PH_UP)
    print(f"슬라이더 최고 높이 {q[1]:.3f} m, 지렛대 코(축에서 270 mm) 윗변 최고 높이 {s.O[1]+0.27*math.sin(s.PH_UP)+s.HW:.3f} m, 보조 스프링 최소 길이 {s.Ls(s.U0+0.105,0)*1000:.0f} mm")
    print(f"== 바닥 발: 피벗 {PF}, 길이 {LF*1000:.0f} mm, 발 위치 z = {ZFOOT:.3f} (뒷바퀴 접지점에서 {abs(ZFOOT-ZR)*1000:.0f} mm 뒤), 필요한 바닥 마찰 tan = {math.tan(LAM):.2f}")
    ph_td = bis(lambda ph: (leg_angle(ph) > LAM + 1e-9) - 0.5, 0.0, s.PH_UP)
    print(f"  들어올림: 크랭크 {RC*1000:.0f} mm, 귀 {RL*1000:.0f} mm, 막대 {LROD*1000:.1f} mm, 장공 {(LROD-dist(crank_pin(s.PH_UP),lug(LAM)))*1000:.1f} mm 이상. 좌판 수평에서 다리 {math.degrees(leg_angle(0)):.0f} deg (발 높이 {(PF[1]-LF*math.cos(leg_angle(0)))*1000:.0f} mm, 발 끝 z = {PF[0]-LF*math.sin(leg_angle(0)):.3f})")
    print(f"  발이 바닥에 닿는 때: 지렛대 {math.degrees(ph_td):.1f} deg = 좌판 {math.degrees(s.th_of(ph_td)):.1f} deg | 지렛대 2/4/6 deg 에서 다리 각 " + ' / '.join(f"{math.degrees(leg_angle(math.radians(a))):.1f}" for a in (2,4,6)))
    for mus in (0.2, 0.3, 0.4):
        row = []
        for m in (40, 70, 120):
            f33 = foot_state(s.TH_UP, m, mus)
            tip = limit(lambda th: foot_state(th, m, mus)['front'] >= 0); vault = limit(lambda th: foot_state(th, m, mus)['rear'] >= 0)
            row.append(f"{m} kg: H={f33['H']:.0f} N, 다리 수직력 {f33['Fv']:.0f} N, 뒤로 들림 한계 {tip:.1f}, 뒷바퀴 뜸 한계 {vault:.1f}")
        print(f"좌판 마찰 {mus}: " + ' | '.join(row))
