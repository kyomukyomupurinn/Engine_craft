import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

# --- 1. カスタム燃料の登録 ---
paraffin_card = """
fuel Paraffin_Koudai  C 32   H 66   wt%=100.0
h,cal=-1800.0      t,k=298.15   rho=0.93
"""
add_new_fuel('Paraffin_Koudai', paraffin_card)
cea = CEA_Obj(oxName='LOX', fuelName='Paraffin_Koudai', pressure_units='MPa', cstar_units='m/sec', temperature_units='K')

# --- 2. パラメータ設定 ---
Pc = 3.0
eps = 10.0
# O/Fを 1.0 から 8.0 まで 50分割した配列を作成
MR_array = np.linspace(1.0, 8.0, 50)

isp_list = []
cstar_list = []

# --- 3. ループ計算 ---
print("計算中...")
for MR in MR_array:
    isp_sl, mode = cea.estimate_Ambient_Isp(Pc=Pc, MR=MR, eps=eps, Pamb=0.101325)
    cstar = cea.get_Cstar(Pc=Pc, MR=MR)
    isp_list.append(isp_sl)
    cstar_list.append(cstar)

# --- 4. グラフ描画 ---
fig, ax1 = plt.subplots(figsize=(8, 5))

# 比推力のプロット (左軸)
ax1.plot(MR_array, isp_list, 'b-', label='Isp (sec)')
ax1.set_xlabel('Mixture Ratio (O/F)')
ax1.set_ylabel('Specific Impulse Isp (sec)', color='b')
ax1.tick_params('y', colors='b')
ax1.grid(True)

# 特性排気速度のプロット (右軸)
ax2 = ax1.twinx()
ax2.plot(MR_array, cstar_list, 'r--', label='c* (m/s)')
ax2.set_ylabel('Characteristic Velocity c* (m/s)', color='r')
ax2.tick_params('y', colors='r')

plt.title('LOX / Paraffin Performance vs Mixture Ratio')
fig.tight_layout()
plt.show() # グラフのウィンドウを表示