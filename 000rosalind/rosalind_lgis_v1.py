from itertools import combinations
#operator.le()

x = "5 1 4 2 3".split()
perms = []
y = []

for i in range(2, len(x)+1):
    for c in combinations(x, i):
        perms.append("".join(c))
for i in range(0, len(perms)):   
    x = list(perms[i])    
    if (x == sorted(x)):
        y.append("".join(str(x)))
print(max((y), key=len))
y = []
for i in range(0, len(perms)):   
    x = list(perms[i])
    if (x == sorted(x)[::-1]):
        y.append("".join(str(x)))
print(max((y), key=len))