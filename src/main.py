from utils import *

def Q1Full(time,N,X0,X1,T0,iterations,question):

    experimental = Q1_Q2(time,N,X0,X1,T0)
    proof = proof_Q1(time,N,X0,X1,T0,iterations)
    graph(time,proof,experimental,N,question)


def Q2Full(time,N,L,X0,X1,T0,iterations,question):

    experimental = Q1_Q2(time,N,X0,X1,T0)
    proof = proof_Q2(time,N,X0,X1,T0,iterations)
    graph(time,proof,experimental,N,question)


def Q3Full(time,N,L,X0,X1,question):

    experimental = Q3(time,N,L,X0,X1)
    proof = proof_Q3(time,N,L,X0,X1)
    graph(time,proof,experimental,L,N,question)


if __name__ == "__main__":
    Q1Full(time=10,N=15,X0=0,X1=0,T0=1,iterations=30,question="Q1")
    print('Question 1 was just plotted!')

    Q2Full(time=10,N=15,X0=1,X1=0,T0=0,iterations=30,question="Q2")
    print('Question 2 was just plotted!')

    Q3Full(time=10,N=25,X0=0,X1=0,question="Q3")
    print('Question 3 was just plotted!')