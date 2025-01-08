import numpy as np

omegaRF = 325e6
periodRF = 1/omegaRF
print(periodRF)

steps = 5

vals = np.linspace(0,periodRF,num=steps)

for i in range(steps):
    print("param beamtime="+str(vals[i]))

''' 20 steps:
param beamtime=0.0
param beamtime=1.619433198380567e-10
param beamtime=3.238866396761134e-10
param beamtime=4.8582995951417e-10
param beamtime=6.477732793522268e-10
param beamtime=8.097165991902835e-10
param beamtime=9.7165991902834e-10
param beamtime=1.1336032388663968e-09
param beamtime=1.2955465587044535e-09
param beamtime=1.4574898785425102e-09
param beamtime=1.619433198380567e-09
param beamtime=1.7813765182186237e-09
param beamtime=1.94331983805668e-09
param beamtime=2.105263157894737e-09
param beamtime=2.2672064777327936e-09
param beamtime=2.4291497975708505e-09
param beamtime=2.591093117408907e-09
param beamtime=2.7530364372469635e-09
param beamtime=2.9149797570850204e-09
param beamtime=3.076923076923077e-09
'''

''' 5 steps:
param beamtime=0.0
param beamtime=7.692307692307692e-10
param beamtime=1.5384615384615385e-09
param beamtime=2.3076923076923076e-09
param beamtime=3.076923076923077e-09
'''