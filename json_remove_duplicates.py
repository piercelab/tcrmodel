#!/usr/bin/env python
from sys import argv
import json
script, f1, f2 = argv

jsonFile = open(f1, 'r')
data = json.load(jsonFile)

size=len(data)
print size

values = [];
uniqueNames = [];
for i in data:
    if(i["name"] not in uniqueNames):
         uniqueNames.append(i["name"]);
         values.append(i)
jsonFile.close()


print len(values)
with open(f2, 'w') as outfile:
    json.dump(values, outfile)
