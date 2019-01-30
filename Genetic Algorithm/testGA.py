from numpy import *
initialPop=6
pc = .25
geneNum= 4
mutationRate = .1; #pm
def eval(chromosomes):
    fobj=[]
    for i in range(0,len(chromosomes)):
        [a,b,c,d]= chromosomes[i]
        func = abs((a + 2*b+3*c+4*d)-30)
        fobj.append(func)
    return fobj

def selection(fobj):
    fitness=[1/(1.0+fobj[i]) for i in range(0,size(fobj))]
    return fitness
def probability(fitness):
    total = sum(fitness)
    p=[fitness[i]/total for i in range(0,size(fitness))]
    return p
def cumProb(P):
    C=[sum(P[0:i]) for i in range(1,size(P)+1)]
    return C
def nextPop(chromosomes,C,R):
    newChromosomes = chromosomes[0]
    for i in range(0,size(R)):
        res=where(R[i]>array(C))
        res=len(res[0])
        if i==0:
            newChromosomes = chromosomes[res]
        else:
            newChromosomes=vstack((newChromosomes,chromosomes[res]))
    return newChromosomes

def crossover(newChromosomes):
    k=0
    parent=[]
    parentNums=[]
    R=random.rand(initialPop)   
    while k<len(newChromosomes):
        if R[k]<pc:
            parent.append(newChromosomes[k])
            parentNums.append(k)
        k=k+1
    C=random.randint(1,geneNum-1,len(parent))
    for i in range(0,len(parentNums)):
        crossedParent =[parent[i][0:C[i]],parent[(i+1)%(len(parent))][C[i]:geneNum]]
       #flatten
        flatParent =[val for sublist in crossedParent for val in sublist]
        newChromosomes[parentNums[i]]=flatParent
    return newChromosomes   

def mutation(Chromosomes):
    totGeneLength = initialPop*geneNum
    mutNum=random.rand(totGeneLength)
    muts=where(mutNum<mutationRate)
    newNums = random.randint(31,size=len(muts[0]))
 
    if(len(muts[0])>0):
        for i in range(0,len(muts[0])):
            ind1=int(muts[0][i]/(geneNum));
            if ind1<0:
                ind1=0;
            if muts[0][i]%geneNum==0:
                ind1=ind1-1
            ind2=muts[0][i]-ind1*geneNum-1
            # print("ind",i,ind1,ind2)
            Chromosomes[ind1][ind2]=newNums[i]
            # print("New Chrom",Chromosomes[ind1])
    return Chromosomes
    
def main():
#create initial chromosomes
    initialPop=6
    cr = .25
    chromosomes=random.randint(31,size=(initialPop,geneNum))
    print("initial",chromosomes)

    rangeMax = 100
    for i in range(0,rangeMax): 
        fobj=eval(chromosomes)
        fitness=selection(fobj)
        P=probability(fitness)
        C=cumProb(P)
        R=random.rand(initialPop)
        newPop = nextPop(chromosomes,C,R)
        crossOverPop=crossover(newPop)
        mutatedPop=mutation(crossOverPop)
        chromosomes = mutatedPop
        bestChromosome = fitness.index(min(fitness))
        bc =chromosomes[bestChromosome]
        result= bc[0]+2*bc[1]+3*bc[2]+4*bc[3]
        print("Iteration %d best chromosome: %s\nvalue: %d\n")%(i,bc,result)
        
        
  
    bestChromosome = fitness.index(min(fitness))
    bc = chromosomes[bestChromosome]
    result= bc[0]+2*bc[1]+3*bc[2]+4*bc[3]
    
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\nFinal Result After %d iterations: %s\nvalue: %d\n")%(rangeMax,bc,result)
main()