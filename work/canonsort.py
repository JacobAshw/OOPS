import random
import copy
letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']
lettersHat = ['a^', 'b^', 'c^', 'd^', 'e^', 'f^', 'g^', 'h^', 'i^', 'j^', 'k^']

lettersCap = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

def gen_random_pot(leng, commas):
    pot = ""
    for i in range(leng):
        roll = random.choice([0,1])
        if(roll):
            pot = pot + random.choice(lettersCap)
        else:
            pot = pot + random.choice(letters)
    commalocs = []
    for i in range(commas):
        newloc = random.randint(1, leng-1)
        while(newloc in commalocs):
            newloc = random.randint(1, leng-1)
        commalocs.append(newloc)
    commalocs.sort()
    commalocs.append(leng)
    finalpot = pot[0:commalocs[0]]
    for i in range(len(commalocs)-1):
        finalpot = finalpot + "," + pot[commalocs[i]:commalocs[i+1]]

    return finalpot
    

def str_to_pot(str: str):
    pottmp = str.split(",")
    pot = []
    for tile in pottmp:
        tl = []
        for char in tile:
            chara = char
            if(chara not in letters):
                chara = lettersHat[letters.index(char.lower())]
            tl.append(chara)
        pot.append(tl)
    return pot


def hat(char: str):
    if(char in letters):
        return lettersHat[letters.index(char)]
    return letters[lettersHat.index(char)]

def shuffle_pot(pot):
    ltrs = []
    for tile in pot:
        for char in tile:
            if(char[0] not in ltrs):
                ltrs.append(char[0])
    ltrs2 = copy.deepcopy(ltrs)
    random.shuffle(ltrs2)
    hatswap = [random.randint(0,1) for l in ltrs] 

    swap = {}
    for i in range(len(ltrs)):
        if(hatswap[i] == 1):
            swap.update({ltrs[i] : hat(ltrs2[i])})
            swap.update({hat(ltrs[i]) : ltrs2[i]})
        else:
            swap.update({ltrs[i] : ltrs2[i]})
            swap.update({hat(ltrs[i]) : hat(ltrs2[i])})

    newpot = []
    for tile in pot:
        newtile = []
        for char in tile:
            newtile.append(swap.get(char))
        newpot.append(newtile)
    random.shuffle(newpot)
    return newpot


def alpha(pot):
    letterseen = 0
    alphamap = {}
    newpot = []
    for tile in pot:
        newtile = []
        for char in tile:
            if(char not in alphamap):
                alphamap.update({char : letters[letterseen]})
                alphamap.update({hat(char) : hat(letters[letterseen])})
                letterseen = letterseen + 1
            newtile.append(alphamap.get(char))
        newtile.sort()
        newpot.append(newtile)
    return newpot

def cannonsort(pot):
    pot = sorted(pot, key=lambda a:(len(a)))
    # if(1==1):
    #     return pot
    # for tile in pot:
    #     tile.sort(reverse=True)
    pot = alpha(pot)
    l = 1
    for tile in pot:
        l = l + len(tile)

    letterhigh = {}
    for tile in pot:
        for chara in tile:
            char = chara[0]
            if(char not in letterhigh):
                letterhigh.update({char : 0})
            letterhigh.update({char : letterhigh.get(char) + pow(l, len(tile))})
    hathigh = {}
    for tile in pot:
        for char in tile:
            if(char not in hathigh):
                hathigh.update({char : 0})
            hathigh.update({char : hathigh.get(char) + pow(l, len(tile))})
    highs = sorted(hathigh.keys(), key=lambda a:(letterhigh.get(a[0]), hathigh.get(a)), reverse=True)

    lengths = [len(tile) for tile in pot]
    lengths = list(set(lengths))
    lengths.sort()
    oflen = {}
    for length in lengths:
        oflen.update({length : []})
    
    newpot = []
    for tile in pot:
        newtile = sorted(tile, key=lambda a:highs.index(a))
        oflen.get(len(newtile)).append(newtile)

    lengths = reversed(lengths)
    for lent in lengths:
        tiles = oflen.get(lent)
        newtiles = sorted(tiles, key=lambda a:([highs.index(b) for b in a]))
        newpot.extend(newtiles)


    # newpot = sorted(newpot, key=lambda a:(100000-len(a), a))
    return newpot

sizemin = 8
sizemax = 20

commasmin = 2
commasmax = 6

successes = 0
tries = 0

for iter in range(1000000):
    size = random.randint(sizemin, sizemax)
    commas = random.randint(commasmin, commasmax)
    pot = str_to_pot(gen_random_pot(size, commas))
    
    others = []
    for jter in range(1):
        others.append(shuffle_pot(pot))
    
    cannon = cannonsort(pot)
    othercannons = [cannonsort(ot) for ot in others]
    
    match = True
    for i in range(1):
        match = match and (str(alpha(cannon)) == str(alpha(othercannons[i])))
    # print(match)

    tries = tries + 1
    if(not match):
        print("Miss-Match")
        print(pot)
        print(cannon)
        print(alpha(cannon))
        for other in othercannons:
            print(others[othercannons.index(other)])
            print(other)
            print(alpha(other))
    else:
        successes = successes + 1

    # if(successes % 100 == 0):
    print("Cannon worked " + str(successes) + " of " + str(tries) + " tries (" + str(successes/tries*100) + "%)")
    # print("------------------------------")



# pot = str_to_pot("cAb,bA,aCa,cC")
# pot = str_to_pot("cd,C,D")

# others = []
# for jter in range(10):
#     others.append(shuffle_pot(pot))

# cannon = cannonsort(pot)
# othercannons = [cannonsort(ot) for ot in others]

# match = True
# for i in range(10):
#     match = match and (str(cannon) == str(othercannons[i]))
# # print(match)
# if(not match):
#     print("FUCK")
#     print(pot)
#     print(cannon)
#     for other in othercannons:
#         print(others[othercannons.index(other)])
#         print(other)