import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

# --- 1. ABS樹脂（カスタム燃料）の定義と登録 ---
# スクリプト実行ごとに毎回この定義を読み込ませる必要があります
abs_card = """
fuel ABS_HyperTEK  C 6.67   H 8.41   N 0.81   wt%=100.0
h,cal=59395.0      t,k=298.15   rho=1.04
"""
add_new_fuel('ABS_HyperTEK', abs_card)

# --- 2. RocketCEA設定 ---
cea = CEA_Obj(oxName='N2O', fuelName='ABS_HyperTEK', pressure_units='MPa', cstar_units='m/sec', temperature_units='K')

# --- 3. 改良版エンジン諸元 ---
m_ox_dot = 0.117          # 酸化剤流量 (kg/s)
D_p_init = 0.025          # 初期ポート径 (m)
L_p = 0.30                # ポート長さ (m)
rho_f = 1040.0            # ABS密度
D_t = 0.010               # ノズルスロート径 10 mm
eps = 3.0                 # 最適化したノズル膨張比 (出口径 20 mm相当)

# スロート断面積
A_t = (np.pi * D_t**2) / 4.0

# 後退速度定数
a = 4.4e-5 
n = 0.5
burn_time = 2.4
dt = 0.01
time_array = np.arange(0, burn_time + dt, dt)

thrust_list, pc_list, of_list = [], [], []
D_p = D_p_init

print("シミュレーション実行中（スロート10mm仕様）...")
for t in time_array:
    # 燃料質量流量の計算
    A_p = (np.pi * D_p**2) / 4.0
    G_ox = m_ox_dot / A_p
    r_dot = a * (G_ox ** n)
    m_f_dot = rho_f * np.pi * D_p * L_p * r_dot
    m_total_dot = m_ox_dot + m_f_dot
    O_F = m_ox_dot / m_f_dot
    
    # 燃焼室圧力 Pc の動的計算 (仮のc*を1400 m/sとして逆算)
    cstar_est = 1400.0 
    Pc = (m_total_dot * cstar_est) / A_t / 1e6  # MPaに変換
    
    # CEAによる性能評価
    isp_sl, mode = cea.estimate_Ambient_Isp(Pc=Pc, MR=O_F, eps=eps, Pamb=0.101325)
    
    # 推力の計算 (効率85%を考慮)
    thrust = (isp_sl * 0.85) * m_total_dot * 9.80665
    
    thrust_list.append(thrust)
    pc_list.append(Pc)
    of_list.append(O_F)
    D_p += 2.0 * r_dot * dt

# グラフ描画
fig, ax1 = plt.subplots(figsize=(8, 5))
ax1.plot(time_array, thrust_list, 'b-', linewidth=2, label='Thrust (N)')
ax1.set_ylim(0, 500)
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('Thrust (N)', color='b')
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(time_array, pc_list, 'g--', linewidth=2, label='Chamber Pressure (MPa)')
ax2.set_ylabel('Chamber Pressure (MPa)', color='g')
ax2.set_ylim(0, 3.5)

plt.title('Thrust & Pressure History (Throat: 10mm, eps: 4.0)')
fig.tight_layout()
plt.show()

print(f"--- 10mmスロート シミュレーション結果 ---")
print(f"平均推力: {np.mean(thrust_list):.1f} N")
print(f"平均燃焼室圧力: {np.mean(pc_list):.2f} MPa")