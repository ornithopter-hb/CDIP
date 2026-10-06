"""좌판 최대 올림각 검토 - 준정적 계산 (표준 라이브러리만 사용)
(1) 최소각: 좌판이 theta 만큼 올라간 자세에서, 몸통 경사 제한 안에서 질량중심(CoM)을 발 위에 올려 좌판을 떠날 수 있는가
    인체 모델은 기립역학_모델.py 와 같은 분절 길이·질량비 (신장 1.70 m, 체중 70 kg). 엉덩이는 좌판에 붙어 앞 피벗을 중심으로 같이 돈다고 봄.
(2) 최대각: 기울어진 좌판에 앉은 사용자가 휠체어를 뒤로 미는 힘 H 때문에 휠체어가 미끄러지거나 앞바퀴가 들리지 않는가 (브레이크를 건 상태)
좌표: x = 전방(+), y = 위(+), 단위 m, N.  실행: python3 좌판각도_검토_계산.py
"""
import math
g=9.81; H=1.70; M=70.0
h_a=0.039*H; L_s=(0.285-0.039)*H; L_t=(0.530-0.285)*H; L_T=(0.818-0.530)*H
L_heel=0.0454/1.74*H; L_toe=0.2018/1.74*H
m_f,m_s,m_th,m_T=2*0.0145*M,2*0.0465*M,2*0.100*M,0.678*M
c_s,c_th,c_T=0.433,0.433,0.626
P=(0.19,0.455)
def bis(f,lo,hi):
    flo=f(lo)
    for _ in range(60):
        m=0.5*(lo+hi); fm=f(m)
        if (flo>0)==(fm>0): lo,flo=m,fm
        else: hi=m
    return 0.5*(lo+hi)
class Body:
    def __init__(s,bhip,a0deg):
        s.hip0=(P[0]-bhip,0.47+0.09); a0=math.radians(a0deg)
        cb=(s.hip0[1]-h_a-L_s*math.cos(a0))/L_t; b0=-math.acos(cb)
        s.A=(s.hip0[0]-L_s*math.sin(a0)-L_t*math.sin(b0),h_a)
        s.heel=s.A[0]-L_heel; s.toe=s.A[0]+L_toe
        s.Cf=(s.heel+0.5*(L_heel+L_toe),0.5*h_a)
    def posture(s,th):
        x=s.hip0[0]-P[0]; y=s.hip0[1]-P[1]
        Hp=(P[0]+x*math.cos(th)+y*math.sin(th), P[1]-x*math.sin(th)+y*math.cos(th))
        dx=Hp[0]-s.A[0]; dy=Hp[1]-s.A[1]; d=math.hypot(dx,dy)
        if d>L_s+L_t: return None
        al=math.atan2(dx,dy); be=math.acos((L_s**2+d**2-L_t**2)/(2*L_s*d)); a=al+be
        K=(s.A[0]+L_s*math.sin(a),s.A[1]+L_s*math.cos(a))
        return Hp,K,a
    def stat(s,th,target='heel',gm_fix=None):
        r=s.posture(th)
        if r is None: return None
        Hp,K,a=r
        Cs=(K[0]+c_s*(s.A[0]-K[0]),K[1]+c_s*(s.A[1]-K[1]))
        Cth=(Hp[0]+c_th*(K[0]-Hp[0]),Hp[1]+c_th*(K[1]-Hp[1]))
        def com(gm):
            CT=(Hp[0]+c_T*L_T*math.sin(gm),Hp[1]+c_T*L_T*math.cos(gm))
            return ((m_f*s.Cf[0]+m_s*Cs[0]+m_th*Cth[0]+m_T*CT[0])/M,(m_f*s.Cf[1]+m_s*Cs[1]+m_th*Cth[1]+m_T*CT[1])/M), CT
        tx={'heel':s.heel,'ankle':s.A[0]}[target]
        if gm_fix is not None: gm=gm_fix
        elif com(math.radians(89))[0][0]<tx: return dict(gm=None)
        elif com(math.radians(-30))[0][0]>tx: gm=math.radians(-30)
        else: gm=bis(lambda q: com(q)[0][0]-tx, math.radians(-30), math.radians(89))
        c,CT=com(gm)
        Mk=g*(m_T*(K[0]-CT[0])+m_th*(K[0]-Cth[0]))/M
        Mh=g*m_T*(CT[0]-Hp[0])/M
        kflex=math.degrees(a-math.atan2(Hp[0]-K[0],Hp[1]-K[1]))
        return dict(gm=gm,Mk=Mk,Mh=Mh,hipy=Hp[1],hipx=Hp[0],kflex=kflex,a=math.degrees(a),comy=c[1],comx=c[0],K=K,Hp=Hp)
if __name__=='__main__':
    # ---- (1) 최소각: 몸통 경사 제한 아래에서 준정적으로 좌판을 떠날 수 있는가 ----
    print('== 최소각 ==')
    for bhip,a0 in ((0.25,25),(0.25,20),(0.25,15),(0.13,25)):
        B=Body(bhip,a0)
        def lean(d,target): 
            r=B.stat(math.radians(d),target); 
            return None if (r is None or r.get('gm') is None) else math.degrees(r['gm'])
        out=[]
        for target in ('heel','ankle'):
            for lim in (30.0,20.0):
                th=None
                for d10 in range(0,701):
                    v=lean(d10/10,target)
                    if v is not None and v<=lim: th=d10/10; break
                out.append((target,lim,th))
        r20=B.stat(math.radians(20),gm_fix=math.radians(30)); r0=B.stat(0,gm_fix=math.radians(30))
        print(f'고관절 {bhip*100:.0f} cm 뒤, 배측굴곡 {a0}:',' | '.join(f'CoM {t} 위·몸통 {l:.0f} 이하 -> {th}' for t,l,th in out))
        for d in (0,20,30,40,50):
            r=B.stat(math.radians(d),gm_fix=math.radians(30)); ra=B.stat(math.radians(d),'ankle')
            print(f"   {d}: 몸통 30 고정 -> CoM 뒤꿈치 기준 {(r['comx']-B.heel)*100:+.1f} cm, 무릎 {r['Mk']:.2f} ({r['Mk']/2.07*100:.0f} % of 2.07), 고관절 {r['Mh']:.2f}, 배측굴곡 {r['a']:.1f}, 고관절 높이 +{(r['hipy']-0.56)*100:.1f}; CoM 발목 위에 필요한 몸통 {('%.1f'%math.degrees(ra['gm'])) if ra and ra.get('gm') is not None else '없음'}")
    # ---- (2) 최대각: 휠체어가 밀리거나 뒤로 들리지 않는가 ----
    print('== 최대각 ==')
    zr,zf=-0.17,0.23      # 뒷바퀴·앞바퀴 접지점 (휠체어 3D 모델 값, 앞바퀴 위치는 임의값)
    Wc=20*g; zc=0.0        # 휠체어 무게 20 kg (사양표), 무게중심 위치는 가정
    b=0.25                 # 앞 피벗-하중점 (v3)
    def chair(th,m,mus):
        V=0.7*m*g; Hh=max(0.0,V*math.tan(th-math.atan(mus)))
        zs=P[0]-b*math.cos(th); hs=0.47+b*math.sin(th)
        Rr=(V*(zf-zs)+Wc*(zf-zc)+Hh*hs)/(zf-zr); Rf=V+Wc-Rr
        return Hh,Rr,Rf
    def lim(f):
        for d10 in range(0,900):
            if not f(math.radians(d10/10)): return (d10-1)/10
        return 90.0
    for mus in (0.0,0.2,0.3,0.4,0.5):
        row=[]
        for m in (40,70,120):
            tip=lim(lambda th: chair(th,m,mus)[2]>=0)
            s3=lim(lambda th: chair(th,m,mus)[0]<=0.3*chair(th,m,mus)[1])
            s5=lim(lambda th: chair(th,m,mus)[0]<=0.5*chair(th,m,mus)[1])
            s7=lim(lambda th: chair(th,m,mus)[0]<=0.7*chair(th,m,mus)[1])
            row.append(f'{m}kg: 들림 {tip:.1f}, 미끄럼(0.3/0.5/0.7) {s3:.1f}/{s5:.1f}/{s7:.1f}')
        print(f'좌판 마찰 {mus}: H=0 한계 {math.degrees(math.atan(mus)):.1f} | '+' | '.join(row))
    for th in (20,30,40):
        H_,Rr,Rf=chair(math.radians(th),70,0.3); print(th,'70kg mu_s 0.3: H',round(H_),'Rr',round(Rr),'Rf',round(Rf),'필요 바퀴 마찰',round(H_/Rr,2))

# ---- (3) 뒤로 들림 한계: 회전 중심을 뒷바퀴 접지점으로 볼 때와 뒷바퀴 축으로 볼 때 ----
def tip_limit(m, mus, pivot_y, zc=0.0, b=0.25, frac=0.7, Wc=20 * g, zr=-0.17):
    """앞바퀴 반력이 0 이 되는 좌판 각 [deg]. pivot_y = 0 (접지점) 또는 0.295 (뒷바퀴 축 높이)
    뒤로 넘기는 모멘트 H*(h_s - pivot_y) 가 버티는 모멘트 V*(z_s - zr) + Wc*(zc - zr) 를 넘는 각"""
    V = frac * m * g
    for d10 in range(0, 900):
        th = math.radians(d10 / 10)
        Hh = max(0.0, V * math.tan(th - math.atan(mus)))
        zs = P[0] - b * math.cos(th); hs = 0.47 + b * math.sin(th)
        if Hh * (hs - pivot_y) > V * (zs - zr) + Wc * (zc - zr): return (d10 - 1) / 10
    return 90.0
if __name__ == '__main__':
    print('== 뒤로 들림 한계 (브레이크 걸림) ==')
    for mus in (0.0, 0.2, 0.3, 0.4, 0.5):
        print(f'좌판 마찰 {mus}: ' + ' | '.join(f"{m} kg: 접지점 {tip_limit(m, mus, 0.0):.1f}, 바퀴축 {tip_limit(m, mus, 0.295):.1f}" for m in (40, 70, 120)))
    print('휠체어 무게중심 위치 민감도 (120 kg, 좌판 마찰 0.3): ' + ' | '.join(f"zc={zc:+.2f}: 접지점 {tip_limit(120, 0.3, 0.0, zc):.1f}, 바퀴축 {tip_limit(120, 0.3, 0.295, zc):.1f}" for zc in (-0.05, 0.0, 0.05)))
    print('하중점 위치 민감도 (120 kg, 좌판 마찰 0.3): ' + ' | '.join(f"b={bb*100:.0f} cm: 접지점 {tip_limit(120, 0.3, 0.0, b=bb):.1f}, 바퀴축 {tip_limit(120, 0.3, 0.295, b=bb):.1f}" for bb in (0.13, 0.20, 0.25, 0.30)))
    for th in (30, 35, 40):
        V = 0.7 * 120 * g; t = math.radians(th); Hh = V * math.tan(t - math.atan(0.3)); zs = P[0] - 0.25 * math.cos(t); hs = 0.47 + 0.25 * math.sin(t)
        print(f"{th} deg, 120 kg, 좌판 마찰 0.3: H={Hh:.0f} N, 넘기는 모멘트 접지점 {Hh*hs:.0f} / 바퀴축 {Hh*(hs-0.295):.0f} N m, 버티는 모멘트 {V*(zs+0.17)+20*g*0.17:.0f} N m")
