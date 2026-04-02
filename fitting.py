import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

# --- 1. ABS樹脂の定義とCEA準備 ---
abs_card = """
fuel ABS_HyperTEK  C 6.67   H 8.41   N 0.81   wt%=100.0
h,cal=59395.0      t,k=298.15   rho=1.04
"""
try:
    add_new_fuel('ABS_HyperTEK', abs_card)
except Exception:
    pass 
cea = CEA_Obj(oxName='N2O', fuelName='ABS_HyperTEK', pressure_units='MPa', cstar_units='m/sec', temperature_units='K')

# --- 2. 実測データの読み込み ---
df = pd.read_csv('thrust_COLOURS.csv')
time_meas = df['Time'].values
thrust_meas = df['Thrust'].values

# 点火タイミングの同期
threshold = 10.0 
start_idx = np.where(thrust_meas > threshold)[0][0]
time_meas = time_meas[start_idx:] - time_meas[start_idx]
thrust_meas = thrust_meas[start_idx:]
burn_time = 2.73 

# --- 3. チューニングパラメータ ---
D_p = 0.0093         # 初期ポート径 (m)　   実測データの形状に合わせて調整
L_p = 0.240         # グレイン長さ (m)　   実測データの形状に合わせて調整
rho_f = 1040.0      # ABS密度 (kg/m^3)　   ABSの素材特性に基づく
m_ox_dot = 0.122     # 平均酸化剤流量 (kg/s) ※実測データから逆算して調整0.334/2.73
Pc = 2.0           # 平均燃焼室圧力 (MPa) これは勘
eps = 3.24          # ノズル膨張比  実測データの形状に合わせて調整
eta_total = 0.869    # 総合効率

a = 0.769e-4        # 定数 a (真値)     # 定数 a
n = 0.45            # べき乗 n

# --- 4. 内部弾道シミュレーションループ ---
dt = 0.01
time_sim = np.arange(0, burn_time, dt)
thrust_sim = []
O_F_list = []
D_p_list = []
Tc_list = []
m_f_total = 0.0

g0 = 9.80665

for t in time_sim:
    # 幾何学計算
    A_p = (np.pi * D_p**2) / 4.0
    G_ox = m_ox_dot / A_p
    
    # 後退速度と燃料流量
    r_dot = a * (G_ox ** n)     #あくまでも経験式に基づく後退速度の計算で平均的な挙動を模倣
    m_f_dot = rho_f * np.pi * D_p * L_p * r_dot #グレインの体積変化から燃料流量を計算
    
    # 熱力学計算 (CEA)
    O_F = m_ox_dot / m_f_dot    #燃料流量からO/Fを計算
    isp_sl, mode = cea.estimate_Ambient_Isp(Pc=Pc, MR=O_F, eps=eps, Pamb=0.101325)
    
    # 推力計算
    thrust = eta_total * isp_sl * (m_ox_dot + m_f_dot) * g0
    
    # 記録と更新
    thrust_sim.append(thrust)
    O_F_list.append(O_F)
    D_p_list.append(D_p * 1000) 
    m_f_total += m_f_dot * dt
    
    D_p += 2.0 * r_dot * dt

# --- 5. 結果の描画 ---
print(f"--- フィッティング結果の評価 ---")
print(f"シミュレーション燃料消費量: {m_f_total * 1000:.1f} g (目標: 46g)")
print(f"最終ポート径: {D_p_list[-1]:.1f} mm")
print(f"O/Fシフト範囲: {O_F_list[0]:.1f} -> {O_F_list[-1]:.1f}")

fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.plot(time_meas, thrust_meas, 'b-', alpha=0.6, linewidth=2, label='Measured Thrust')
ax1.plot(time_sim, thrust_sim, 'r--', linewidth=2, label='Simulated Thrust')
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('Thrust (N)')
ax1.grid(True)
ax1.legend(loc='upper right')

ax2 = ax1.twinx()
ax2.plot(time_sim, O_F_list, 'g:', linewidth=2, label='O/F Ratio (Sim)')
ax2.set_ylabel('O/F Ratio', color='g')
ax2.tick_params('y', colors='g')
ax2.legend(loc='upper left')

plt.title('Thrust Curve Fitting & Internal Ballistics')
fig.tight_layout()
plt.show()

# --- 6. 燃焼ガス物性の代表値計算 ---
# 代表的なO/Fとして，ループで計算されたO/Fの平均値を使用します
MR_avg = np.mean(O_F_list)

# Tc（燃焼室温度）と cstar（特性排気速度）を計算
Tc = cea.get_Tcomb(Pc=Pc, MR=MR_avg)
cstar = cea.get_Cstar(Pc=Pc, MR=MR_avg)
isp_sl, mode = cea.estimate_Ambient_Isp(Pc=Pc, MR=MR_avg, eps=eps, Pamb=0.101325)

# 比熱比と分子量の取得
mw, gamma = cea.get_Chamber_MolWt_gamma(Pc=Pc, MR=MR_avg, eps=eps)

# 気体定数の計算 (R = Ru / M) ※途中計算
# Ru = 8314.46 J/(kmol K)
# R_gas = 8314.46 / mw
Ru = 8314.46
R_gas = Ru / mw

print("\n--- 燃焼性能計算結果 (平均O/Fに基づく) ---")
print(f"燃焼室温度 (Tc) : {Tc:.2f} K")
print(f"特性排気速度 (c*): {cstar:.2f} m/s")
print(f"海面比推力 (Isp) : {isp_sl:.2f} sec")
print(f"ノズルの状態 (Mode): {mode}")
print(f"--- 燃焼ガス物性値 ---")
print(f"平均分子量 (M)  : {mw:.2f} g/mol")
print(f"比熱比 (γ)      : {gamma:.4f}")
print(f"気体定数 (R)    : {R_gas:.2f} J/(kg K)")