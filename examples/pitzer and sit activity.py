from math import sqrt
import pytherm.activity.pitzer_new as pz
import pytherm.activity.sit as sit

c_nacl = 1
c_kcl = 1

phase = {
    'Na_+1': c_nacl,
    'Cl_-1': c_nacl + c_kcl,
    'K_+1': c_kcl,
}


am = pz.Pitzer(phase, dict_mode=True)
y = am.get_y(phase)
print(y)

a = {}
for i in phase:
    a[i] = phase[i] * y[i]
print(a)


am = sit.SIT(list(phase.keys()), dict_mode=True)
y = am.get_y(phase)
print(y)

a = {}
for i in phase:
    a[i] = phase[i] * y[i]
print(a)