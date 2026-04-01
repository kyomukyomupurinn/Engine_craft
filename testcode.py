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


Pc = 2.38     # 燃焼室圧力 (MPa)
MR = 8.16    # 混合比 (O/F)
eps = 4.0   # ノズル膨張比

# 性能計算 (変数の数を2つに修正！)
isp_sl, mode = cea.estimate_Ambient_Isp(Pc=Pc, MR=MR, eps=eps, Pamb=0.101325)
cstar = cea.get_Cstar(Pc=Pc, MR=MR)
Tc = cea.get_Tcomb(Pc=Pc, MR=MR)

print("--- 燃焼性能計算結果 ---")
print(f"燃焼室温度 (Tc) : {Tc:.2f} K")
print(f"特性排気速度 (c*): {cstar:.2f} m/s")
print(f"海面比推力 (Isp) : {isp_sl:.2f} sec")
print(f"ノズルの状態 (Mode): {mode}")