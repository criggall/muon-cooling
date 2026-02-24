beamtime = -0.671
beamstart = -700

toffs0 = beamtime - beamstart/275.89 - 0.07
print(f"toffs0 = {toffs0}")

toffs1 = toffs0 - 0.061 - 0.0026 + 0.0115 + 0.006
print(f"toffs1 = {toffs1}")

toffs2 = toffs1 - 0.064 + 0.149 - 0.0052 - 0.0029
print(f"toffs2 = {toffs2}")

toffs3 = toffs2
print(f"toffs3 = {toffs3}")

toffs4 = toffs3
print(f"toffs4 = {toffs4}")

toffs5 = toffs4 - 0.008
print(f"toffs5 = {toffs5}")

toffs6 = toffs5
print(f"toffs6 = {toffs6}")