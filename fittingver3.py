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

threshold = 10.0 
start_idx = np.where(thrust_meas > threshold)[0][0]
time_meas = time_meas[start_idx:] - time_meas[start_idx]
thrust_meas = thrust_meas[start_idx:]
burn_time = 2.73 

# --- 3. チューニングパラメータとノズル形状 ---
D_p = 0.012         # 初期ポート径 (m)
L_p = 0.240          # グレイン長さ (m)
rho_f = 1040.0       # ABS密度 (kg/m^3)
m_ox_dot = 0.122     # 平均酸化剤流量 (kg/s)
eps_initial = 3.24   # 初期ノズル膨張比
eta_total = 0.869    # 総合効率

a = 0.769e-4         # 定数 a
n = 0.45             # べき乗 n

# スロート径 (m) とエロージョンの設定
D_t_initial = 0.0093          
A_t_initial = np.pi * (D_t_initial**2) / 4.0
A_e = eps_initial * A_t_initial  # ノズル出口面積 (燃焼中一定と仮定)

# 【追加】スロートの半径方向エロージョンレート (m/s)
# ※ 0.05e-3 m/s = 0.05 mm/s (直径では 0.1 mm/s の拡大)
E_r = 0.05e-3  

# --- 4. 内部弾道シミュレーションループ ---
dt = 0.01
time_sim = np.arange(0, burn_time, dt)

# 記録用リスト
thrust_sim, O_F_list, D_p_list = [], [], []
Pc_list, cstar_list, gamma_list, Cf_exp_list = [], [], [], []
eps_list = [] # 【追加】膨張比の記録リスト

m_f_total = 0.0
Pc_current = 2.0  
g0 = 9.80665
D_t_current = D_t_initial

for i, t in enumerate(time_sim):
    # 幾何学計算と流量計算
    A_p = (np.pi * D_p**2) / 4.0
    G_ox = m_ox_dot / A_p
    r_dot = a * (G_ox ** n) 
    m_f_dot = rho_f * np.pi * D_p * L_p * r_dot 
    O_F = m_ox_dot / m_f_dot 
    
    # 【追加】スロート径と膨張比の動的更新
    D_t_current = D_t_initial + 2.0 * E_r * t
    A_t_current = np.pi * (D_t_current**2) / 4.0
    eps_current = A_e / A_t_current
    
    # CEA物性値の取得 (eps_current を使用)
    cstar = cea.get_Cstar(Pc=Pc_current, MR=O_F)
    mw, gamma = cea.get_Chamber_MolWt_gamma(Pc=Pc_current, MR=O_F, eps=eps_current)
    
    # 質量流量とc*から真の燃焼室圧力(Pc)を逆算 (A_t_current を使用)
    Pc_pa = (m_ox_dot + m_f_dot) * cstar / A_t_current
    Pc_current = Pc_pa / 1e6  
    
    # シミュレーション推力の計算 (eps_current を使用)
    isp_sl, mode = cea.estimate_Ambient_Isp(Pc=Pc_current, MR=O_F, eps=eps_current, Pamb=0.101325)
    thrust = eta_total * isp_sl * (m_ox_dot + m_f_dot) * g0
    
    # 実測推力からの推力係数(Cf)の逆算
    if i < len(thrust_meas):
        F_exp = thrust_meas[i]
        Cf_exp = F_exp / (Pc_pa * A_t_current)
    else:
        Cf_exp = np.nan
    
    # 記録
    thrust_sim.append(thrust)
    O_F_list.append(O_F)
    Pc_list.append(Pc_current)
    cstar_list.append(cstar)
    gamma_list.append(gamma)
    Cf_exp_list.append(Cf_exp)
    eps_list.append(eps_current) # 記録
    
    m_f_total += m_f_dot * dt
    D_p += 2.0 * r_dot * dt

# --- 5. 結果の描画 ---
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# グラフ1: 推力とO/F
ax1.plot(time_meas, thrust_meas, 'b-', alpha=0.6, label='Measured Thrust')
ax1.plot(time_sim, thrust_sim, 'r--', label='Simulated Thrust')
ax1.set_ylabel('Thrust (N)')
ax1.grid(True)
ax1.legend(loc='upper left')
ax1_sub = ax1.twinx()
ax1_sub.plot(time_sim, O_F_list, 'g:', linewidth=2, label='O/F Ratio')
ax1_sub.set_ylabel('O/F Ratio', color='g')

# グラフ2: 燃焼室圧力(Pc)と比熱比(γ)
ax2.plot(time_sim, Pc_list, 'k-', linewidth=2, label='Chamber Pressure ($P_c$)')
ax2.set_ylabel('Pressure (MPa)')
ax2.grid(True)
ax2.legend(loc='upper left')
ax2_sub = ax2.twinx()
ax2_sub.plot(time_sim, gamma_list, 'm-.', linewidth=2, label='Specific Heat Ratio ($\gamma$)')
ax2_sub.set_ylabel('Specific Heat Ratio ($\gamma$)', color='m')

# グラフ3: 実測推力係数(Cf)と【追加】ノズル膨張比(ε)
ax3.plot(time_sim, Cf_exp_list, 'c-', linewidth=2, label='Experimental $C_F$')
ax3.set_xlabel('Time (sec)')
ax3.set_ylabel('Thrust Coefficient ($C_F$)')
ax3.grid(True)
ax3.legend(loc='upper left')
ax3_sub = ax3.twinx()
ax3_sub.plot(time_sim, eps_list, 'r:', linewidth=2, label='Expansion Ratio ($\epsilon$)')
ax3_sub.set_ylabel('Expansion Ratio ($\epsilon$)', color='r')
ax3_sub.legend(loc='upper right')

plt.suptitle('Transient Internal Ballistics with Nozzle Erosion')
fig.tight_layout()
plt.show()