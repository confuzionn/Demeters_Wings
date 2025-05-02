###super test zone

import re

strng = "d01 rm 6\n d12 pm 0.0\n d13 ym .0123"
key = "d+[0-9] [rpy]m"
key = "d[0-9]+ [rpy]m [0-9]*.[0-9]*"
x = re.findall(key, strng)
print(x)

# print("{:<14} t".format("asl moment 🖐👇🙌👐👌👉✌"))

#d01 -> d1

x = re.findall("d0[0-9]", strng)
print(x)
for i, heathen in enumerate(x):  
    x[i] = re.sub('0','',heathen)
    strng = re.sub("d0[0-9]", x[i], strng, count=1)
print(strng)

x = re.search("d[0-9]*", strng)
print(x.group())

x = {'yes' : 1}
try: 
    x['yes']
    print("huge w for all of mankind")
except:
    print('NOOOOOOOOOOO')

print(type(-1.0))