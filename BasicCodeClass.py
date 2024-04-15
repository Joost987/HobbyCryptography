import itertools as iter
import math
import random


loweralpha="abcdefghijklmnopqrstuvwxyz"
upperalpha=loweralpha.upper()
totalalpha=loweralpha+upperalpha

f=open("wordlist.txt",'r')
wordlist=f.read().splitlines()

f=open("WarAndPeace.txt",'r')
WarAndPeace=f.read().splitlines()


class Code:
    def __init__(self,encoded,disregardCapitals=True,plain=None):

        if disregardCapitals:
            self.code=encoded.lower()
        else:
            self.code=encoded
        self.codeNoSpaces="".join([char if char!=" " else "" for char in self.code])
        self.len_with_spaces=len(encoded)
        self.len=len(self.codeNoSpaces)
        self.plain=plain
        self.frequencies=None
        

    def FrequencyAnalysis(self):
        
        self.frequencies=dict([(char,0) if chr!=" " else None for char in self.codeNoSpaces])
        for i in range(self.len):
            self.frequencies[self.codeNoSpaces[i]]+=1

    def IndexOfOccurence(self):
        if self.frequencies==None:
            self.FrequencyAnalysis()
        IC=0    
        for letter in self.frequencies.keys():
            IC+=self.frequencies[letter]*(self.frequencies[letter]-1)
        IC=IC/(self.len*(self.len-1))
        return IC

    def NGramsAnalysis(self,N=2,DisregardSpaces=False):
        NGramFrequencies=dict()
        for NGram in iter.product(loweralpha+" ",repeat=N):
            NGramFrequencies.update({"".join(NGram):0})
        
        if DisregardSpaces==True:
            for i in range(self.len-N+1):
                NGramFrequencies.update({self.codeNoSpaces[i:(i+N)]:NGramFrequencies[self.codeNoSpaces[i:i+N]]+1})
        else:
            for i in range(self.len_with_spaces-N+1):
                NGramFrequencies.update({self.code[i:(i+N)]:NGramFrequencies[self.code[i:i+N]]+1})
            
        return NGramFrequencies
    
    def PlainNGramsAnalysis(self,N=2):
        NGramFrequencies=dict()
        for NGram in iter.product(loweralpha+" ",repeat=N):
            NGramFrequencies.update({"".join(NGram):0})
        
        for i in range(len(self.plain)-N+1):
            NGramFrequencies.update({self.plain[i:(i+N)]:NGramFrequencies[self.plain[i:i+N]]+1})
            
        return NGramFrequencies
    
    def PlainnessScore(self,RefTwoGramFreqs,RefThreeGramFreqs): #gives log likeliness
        TwoGramFrequenciesCode=self.PlainNGramsAnalysis()
        ThreeGramFrequenciesCode=self.PlainNGramsAnalysis(N=3)
        sim=0
        a=0.5
        for TwoGram in TwoGramFrequenciesCode.keys():
            #sim+=math.log(RefTwoGramFreqs[TwoGram])*(TwoGramFrequenciesCode[TwoGram])
            sim+=math.log(RefTwoGramFreqs[TwoGram])*(TwoGramFrequenciesCode[TwoGram])
        sim*=a
        for ThreeGram in ThreeGramFrequenciesCode.keys():
            #sim+=math.log(RefTwoGramFreqs[TwoGram])*(TwoGramFrequenciesCode[TwoGram])
            sim+=(1-a)*math.log(RefThreeGramFreqs[ThreeGram])*(ThreeGramFrequenciesCode[ThreeGram])
        
        return sim

        
class RotCypher(Code): #add general alphabet support
    def __init__(self,encoded,alphabet=loweralpha,plain=None,rot=None):
        NotcontainsCaps=set(alphabet).intersection(upperalpha)==set()
        Code.__init__(self,encoded,NotcontainsCaps,plain)
        self.rot=rot
        self.alphabet=alphabet

    def Decode(self):
        if self.rot==None:
            print("Cannot decode without the rotation key, if you want to brute force use DecodeBrute")
            return
        self.plain=""

        for i in range(self.len_with_spaces):
            if self.code[i]==" ":
                self.plain+=" "
            elif self.code[i] in self.alphabet:
                self.plain+=self.alphabet[(self.alphabet.find(self.code[i])-self.rot)%len(self.alphabet)]
            else:
                print("This character \""+self.code[i]+"\" was not in the given alphabet")
                raise(RuntimeError)
    
    def DecodeBrute(self):
        plainlist=[]
        for rot in range(len(self.alphabet)):
            self.rot=rot
            self.Decode()
            plainlist.append(self.plain)
        return plainlist


def EncodeRot(plain,rot,alphabet=loweralpha):
    code=""
    for i in range(len(plain)):
        if plain[i] in alphabet:
            code+=alphabet[(alphabet.find(plain[i])+rot)%len(alphabet)]
        elif plain[i]==" ":
            code+=" "
        else:
            print("Character \""+plain[i]+"\" was not contained in given alphabet")
            raise(RuntimeError)
    return RotCypher(code,alphabet,plain,rot)

class MonoSubCypher(Code):
    def __init__(self,encoded,subKey=None,subAlphabet=None,plainAlphabet=loweralpha,disregardCapitals=False,plain=None):
        Code.__init__(self,encoded,disregardCapitals,plain)
        if subKey!=None and subAlphabet!=None:
            print("You cannot give both a substitution key and a substitution alphabet, give one or the other")
            raise(RuntimeError)
        elif subKey!=None and subAlphabet==None:
            self.subAlpha=""
            for char in subKey:
                if char not in self.subAlpha and char!=" ":
                    self.subAlpha+=char
            for char in plainAlphabet:
                if char not in subKey:
                    self.subAlpha+=char
        else:
            self.subAlpha=subAlphabet
        self.plainAlpha=plainAlphabet
        self.plain=plain
        if self.subAlpha!=None:
            #print(self.subAlpha,len(self.subAlpha),self.plainAlpha,len(self.plainAlpha))
            self.decodeDic=dict([(self.subAlpha[i], self.plainAlpha[i]) for i in range(len(self.subAlpha))])

    def Decode(self):
        if self.decodeDic==None:
            print("Cannot decode without substitution alphabet or substitution key")
            return
        self.plain=""

        for i in range(self.len_with_spaces):
            if self.code[i]==" ":
                self.plain+=" "
            elif self.code[i] in self.subAlpha:
                self.plain+=self.decodeDic[self.code[i]]
            else:
                print("This character \""+self.code[i]+"\" was not in the given alphabet")
                raise(RuntimeError)
        
    def MCMCSolver(self,StartingSubAlpha,RefFreqsTwo,RefFreqsThree,maxiter=50000): #credits to https://maximilianrohde.com/posts/code-breaking-with-metropolis/
        CurrentSubAlpha=StartingSubAlpha
        self.subAlpha=CurrentSubAlpha
        self.decodeDic=dict([(self.subAlpha[i], self.plainAlpha[i]) for i in range(len(self.subAlpha))])
        self.Decode()
        CurrScore=self.PlainnessScore(RefFreqsTwo,RefFreqsThree)

        for i in range(maxiter):

            ProposedSubAlpha=SwapTwoRandom(CurrentSubAlpha)
            self.subAlpha=ProposedSubAlpha
            self.decodeDic=dict([(self.subAlpha[i], self.plainAlpha[i]) for i in range(len(self.subAlpha))])
            self.Decode()
            PropScore=self.PlainnessScore(RefFreqsTwo,RefFreqsThree)

            AcceptProb=min(1,math.exp(PropScore-CurrScore))
            if random.random()<=AcceptProb:
                CurrentSubAlpha=ProposedSubAlpha
                CurrScore=PropScore
            if i%100==0:
                print("Iteration: "+ str(i)+ " " + self.plain)
        print(self.plain)

    def PreMMC(self,RefFreqsTwo,RefFreqsThree,maxiter=50000):
        scoremin=0
        for i in range(30):
            lst=list(loweralpha)
            random.shuffle(lst)
            randAlpha="".join(lst)
            self.subAlpha=randAlpha
            self.decodeDic=dict([(self.subAlpha[i], self.plainAlpha[i]) for i in range(len(self.subAlpha))])
            self.Decode()
            score=self.PlainnessScore(RefFreqsTwo,RefFreqsThree)
            if score<scoremin:
                scoremin=score
                subAlphaMin=self.subAlpha
        return subAlphaMin
            

def EncodeMonoSub(plain,subKey,plainAlphabet=loweralpha):
    subAlpha=dict([(plainAlphabet[i],subKey[i]) for i in range(len(subKey))])
    counter=len(subKey)
    for (i,char) in enumerate(plainAlphabet):
        if char not in subAlpha.values():
            subAlpha.update({plainAlphabet[counter]:char})
            counter+=1
    code=""
    for char in plain:
        if char==" ":
            code+=" "
        else:
            code+=subAlpha[char]
    return MonoSubCypher(code,subKey=subKey,plainAlphabet=plainAlphabet,plain=plain)
    


def PickBestSolution(plainlist,wordlist): #Picks all solutions with most words in the wordlist
    correctnessDict=dict([(plain,sum([1 if word in wordlist else 0 for word in plain.split(" ")])) for plain in plainlist])
    bestScore=max(correctnessDict.values())
    return [plain for plain,score in correctnessDict.items() if score==bestScore]



def PickValidSolution(plainlist,wordlist): #Picks first solution of which all words are in the wordlist
    for plain in plainlist:
        for word in plain.split(" "):
            param=False
            if word not in wordlist:
                break
            param=True
        if param:
            return plain
    return None


def CalculateReferenceFreqs(string,N=2):
    RefFreqs=dict()
    for NGram in iter.product(loweralpha+" ",repeat=N):
        RefFreqs.update({"".join(NGram):1})
    
    string=string.lower()
    filteredstring=""
    for chr in string:
        if chr in (loweralpha+" "):
            filteredstring+=chr
    for i in range(len(filteredstring)-N+1):
        RefFreqs.update({filteredstring[i:i+N]:RefFreqs[filteredstring[i:i+N]]+1})
    
    total=sum(RefFreqs.values())
    for NGram in iter.product(loweralpha+" ",repeat=N):
        RefFreqs.update({"".join(NGram):RefFreqs["".join(NGram)]/total})
    return RefFreqs

def SwapTwoRandom(string):
    indexes=random.sample(range(0,len(string)),k=2)
    lst=list(string)
    lst[indexes[0]],lst[indexes[1]]=lst[indexes[1]],lst[indexes[0]]
    return "".join(lst)



#RefFreqsTwo=CalculateReferenceFreqs("".join(WarAndPeace))
#RefFreqsThree=CalculateReferenceFreqs("".join(WarAndPeace),N=3)