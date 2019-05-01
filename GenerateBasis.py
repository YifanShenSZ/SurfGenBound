# According to user input generate basis.in, specifying the expansion basis functions in SurfGenBound

''' User input '''
DegreeOfFreedom=39
# 1st element in sublist tells the order, others are the degrees following this order
# try 0
    # 2th-order for all
#SelectBasis=[[2,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]]
# try 1
    # 4th-order for carbon skeleton dihedral angle (18,19)
    # 2th-order for others
SelectBasis=[[4,18,19],[2,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]]
# try 2
    # 4th-order for carbon skeleton dihedral angle (18,19) &
    #               C-O bond (1) & C_alpha-H bond (7) & O or H involved C_alpha angle (20-23)
    # 2th-order for others
#SelectBasis=[[4,18,19,1,7,20,21,22,23],[2,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]]
# try 3
    # 4th-order for carbon skeleton bond (2-6) & angle (16,17) & dihedral angle (18,19) &
    #               C-O bond (1) & C_alpha-H bond (7) & O or H involved C_alpha angle (20-23)
    # 2th-order for others
#SelectBasis=[[4,2,3,4,5,6,16,17,18,19,1,7,20,21,22,23],[2,8,9,10,11,12,13,14,15,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]]
# try 4
    # 6th-order for carbon skeleton dihedral angle (18,19) & C-O bond (1) & C_alpha-H bond (7) & O or H involved C_alpha angle (20-23)
    # 4th-order for carbon skeleton bond (2-6) & angle (16,17)
    # 2th-order for others
#SelectBasis=[[6,18,19,1,7,20,21,22,23],[4,2,3,4,5,6,16,17],[2,8,9,10,11,12,13,14,15,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]]

''' Import libraries '''
import sys
import numpy
import scipy.special

''' Auxiliary routines '''
def PickOutOrder(list):
	return list[0]

def AllTerms(order,dof):
    dof.sort()
    temp=numpy.empty((order,int(scipy.special.comb(len(dof)+order-1,order))),dtype=int)
    temp[:,0]=0
    for i in range(1,temp.shape[1]):
        temp[:,i]=temp[:,i-1]
        temp[0,i]=temp[0,i]+1
        for j in range(order-1):
            if(temp[j,i]>len(dof)-1):
                temp[j,i]=0
                temp[j+1,i]=temp[j+1,i]+1
        for j in range(order-2,-1,-1):
            if(temp[j,i]<temp[j+1,i]):
                temp[j,i]=temp[j+1,i]
    output=numpy.empty((order,temp.shape[1]),dtype=int)
    for i in range(output.shape[1]):
        for j in range(order):
            output[j,i]=dof[temp[j,i]]
    return output

''' Do the job '''
temp=0# Check whether all degree of freedoms are assigned correctly
for i in range(len(SelectBasis)):
    temp=temp+len(SelectBasis[i])-1
if(temp!=DegreeOfFreedom):
    print('Program abort: Degree of freedom does not consistent!')
    sys.exit()
SelectBasis.sort(key=PickOutOrder,reverse=True)# Sort orders descendingly
output=[]
for i in range(len(SelectBasis)-1):
    temp=[]
    for j in range(i+1):
        temp=temp+SelectBasis[j][1:]
    temp.sort()
    for j in range(SelectBasis[i][0],SelectBasis[i+1][0],-1):
        output.append(AllTerms(j,temp))
temp=[]
for j in range(len(SelectBasis)):
    temp=temp+SelectBasis[j][1:]
temp.sort()
for j in range(SelectBasis[len(SelectBasis)-1][0],0,-1):
    output.append(AllTerms(j,temp))
with open('basis.in','w') as f:
    order=numpy.empty(len(output),dtype=int)
    order[0]=SelectBasis[0][0]
    for i in range(1,order.shape[0]):
        order[i]=order[i-1]-1
    for i in range(len(output)):
        for j in range(output[i].shape[1]):
            print('%5d'%order[i],file=f,end='')
            for k in range(output[i].shape[0]):
                print('%5d'%output[i][k,j],file=f,end='')
            print(file=f)
    print('%5d'%0,file=f,end='')