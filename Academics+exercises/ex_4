#!/usr/bin/env python
'''
from numpy import loadtxt

lines =loadtxt("ihaveadream.txt", delimiter=" ", unpack=False)

print(lines)
'''

import re

f = open('ihaveadream.txt','r')
x = f.readlines()

wordlist = [re.findall('\w+',x[l]) for l in xrange(len(x))]

'''
for line in wordlist:
    if len(wordlist[line]) == 0:
        wordlist.pop(line)
        line in
'''

#print(wordlist)

word_lis_dict = [wordlist[line][word] for line in xrange(len(wordlist)) for word in xrange(len(wordlist[line]))]

#change every letter to lowercase
word_dict = [re.sub('([A-Z]+)', lambda m: m.group(0).lower(), word_lis_dict[word]) for word in xrange(len(word_lis_dict))]

#print(word_dict)

#diction=dict()

#diction[word_dict[iter]] = [word_dict.count(word_dict[iter]) for iter in xrange(len(word_dict))]

diction=dict()
for iter in xrange(len(word_dict)):
    diction[word_dict[iter]]= word_dict.count(word_dict[iter])
    if word_dict[iter] == 'negro': print(word_dict[iter]+' and the indice '+str(iter))
    
show_lis = sorted(diction, key=diction.__getitem__, reverse=True)

indices = list(xrange(len(show_lis)))

counter = 0
for ind in indices:
    try:
        if len(show_lis[counter]) < 4:
            show_lis.remove(show_lis[counter])
            indices.pop()
        else:
            counter += 1
    except:
        pass

for iter in xrange(10):
    print(show_lis[iter]+': '+str(diction[show_lis[iter]]))

#print(diction['negro'])
#print(wordlist)
#print(wordlist.count('negro'))

'''
count = 0
primary = word_dict[count]
diction = dict()
primary_ind = []
while count != len(word_dict)-1:
    count += 1
    multiplicity = 1
    if len(primary) > 3:
        primary_ind.append(count)
        for word in word_dict:           
            if primary == word:
                diction[primary]=multiplicity+1
                multiplicity += 1
            if word == word_dict[-1]:
                primary = word_dict[count]                
                #print(primary)
    else:
        primary = word_dict[count]
        #print(primary)
#type(diction)
#print(sorted(diction, key=diction.__getitem__, reverse=True))

show_lis = sorted(diction, key=diction.__getitem__, reverse=True)

for iter in xrange(10):
    print(show_lis[iter]+': '+str(diction[show_lis[iter]]))
'''

'''
multiplicity_lis = list(enumerate(diction.values()))
word_lis = list(enumerate(diction.keys()))

print(max(diction.values()))
'''
'''
maxes = []
for iter in xrange(10):
    maxes.append(max(multiplicity_lis))
    
'''



    
'''
wordlist =[re.split('\ ', x[l]) for l in xrange(len(x))]

wordlist =[wordlist[line][word].rstrip('\n') for line in xrange(len(wordlist)) for word in xrange(len(wordlist[line]))]
'''
#print(wordlist)
#print(word)

