from rocketcea.cea_obj_w_units import CEA_Obj

# 酸化剤：液体酸素(LOX)，燃料：ケロシン(RP-1)
cea = CEA_Obj(oxName='LOX', fuelName='RP-1', pressure_units='MPa', cstar_units='m/sec', temperature_units='K')

Pc = 3.0     # 燃焼室圧力 (MPa)
MR = 2.5     # 混合比 (O/F)
eps = 10.0   # ノズル膨張比

# 性能計算 (変数の数を2つに修正！)
isp_sl, mode = cea.estimate_Ambient_Isp(Pc=Pc, MR=MR, eps=eps, Pamb=0.101325)
cstar = cea.get_Cstar(Pc=Pc, MR=MR)
Tc = cea.get_Tcomb(Pc=Pc, MR=MR)

print("--- 燃焼性能計算結果 ---")
print(f"燃焼室温度 (Tc) : {Tc:.2f} K")
print(f"特性排気速度 (c*): {cstar:.2f} m/s")
print(f"海面比推力 (Isp) : {isp_sl:.2f} sec")
print(f"ノズルの状態 (Mode): {mode}")