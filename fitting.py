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
    pass # 既に登録済みの場合はスキップ
cea = CEA_Obj(oxName='N2O', fuelName='ABS_HyperTEK', pressure_units='MPa', cstar_units='m/sec', temperature_units='K')

# --- 2. 実測データの読み込み ---
df = pd.read_csv('thrust_COLOURS.csv')
time_meas = df['Time'].values
thrust_meas = df['Thrust'].values

# 点火タイミングの同期（推力が上がり始めた時間を0秒とする）
threshold = 10.0 # 10Nを超えた時点をスタートとする
start_idx = np.where(thrust_meas > threshold)[0][0]
time_meas = time_meas[start_idx:] - time_meas[start_idx]
thrust_meas = thrust_meas[start_idx:]
burn_time = 2.73 # 燃焼時間 (sec)

# --- 3. ★チューニングパラメータ（ここを弄って波形を合わせる）★ ---
D_p = 0.020         # 初期ポート径 (m) [例: 20mm]
L_p = 0.150         # グレイン長さ (m) [例: 150mm]
rho_f = 1040.0      # ABS密度 (kg/m^3)
m_ox_dot = 0.12     # 平均酸化剤流量 (kg/s) -> 逆算値から推測
Pc = 2.38           # 平均燃焼室圧力 (MPa)
eps = 4.0           # ノズル膨張比
eta_total = 0.80    # 総合効率 -> 逆算値から推測

# 後退速度定数 (r_dot = a * G_ox^n)
a = 0.20e-4         # 定数 a (値を大きくすると推力の傾きが急になる)
n = 0.5             # べき乗 n (通常0.4~0.6)

# --- 4. 内部弾道シミュレーションループ ---
dt = 0.01
time_sim = np.arange(0, burn_time, dt)
thrust_sim = []
O_F_list = []
D_p_list = []
m_f_total = 0.0

g0 = 9.80665

for t in time_sim:
    # 幾何学計算
    A_p = (np.pi * D_p**2) / 4.0
    G_ox = m_ox_dot / A_p
    
    # 後退速度と燃料流量
    r_dot = a * (G_ox ** n)
    m_f_dot = rho_f * np.pi * D_p * L_p * r_dot
    
    # 熱力学計算 (CEA)
    O_F = m_ox_dot / m_f_dot
    isp_sl, mode = cea.estimate_Ambient_Isp(Pc=Pc, MR=O_F, eps=eps, Pamb=0.101325)
    
    # 推力計算
    thrust = eta_total * isp_sl * (m_ox_dot + m_f_dot) * g0
    
    # 記録と更新
    thrust_sim.append(thrust)
    O_F_list.append(O_F)
    D_p_list.append(D_p * 1000) # mm単位で記録
    m_f_total += m_f_dot * dt
    
    D_p += 2.0 * r_dot * dt

# --- 5. 結果の描画 ---
print(f"--- フィッティング結果の評価 ---")
print(f"シミュレーション燃料消費量: {m_f_total * 1000:.1f} g (目標: 40g)")
print(f"最終ポート径: {D_p_list[-1]:.1f} mm")
print(f"O/Fシフト範囲: {O_F_list[0]:.1f} -> {O_F_list[-1]:.1f}")

fig, ax1 = plt.subplots(figsize=(10, 6))

# 実測とシミュレーションの推力比較
ax1.plot(time_meas, thrust_meas, 'b-', alpha=0.6, linewidth=2, label='Measured Thrust')
ax1.plot(time_sim, thrust_sim, 'r--', linewidth=2, label='Simulated Thrust')
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('Thrust (N)')
ax1.grid(True)
ax1.legend(loc='upper right')

# O/Fの推移 (右軸)
ax2 = ax1.twinx()
ax2.plot(time_sim, O_F_list, 'g:', linewidth=2, label='O/F Ratio (Sim)')
ax2.set_ylabel('O/F Ratio', color='g')
ax2.tick_params('y', colors='g')
ax2.legend(loc='upper left')

plt.title('Thrust Curve Fitting & Internal Ballistics')
fig.tight_layout()
plt.show()