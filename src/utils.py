import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from math import log10, floor

################################################################################
#                                                                              #
#                             LOOPED COMPUTING                                 #
#                                                                              #
################################################################################

def Q1_Q2(time,N,L,X0,X1,T0)->np.matrix:
    ''' Returns the temperature of the bar on N intervals of the bar, 
    in different points of time (dt intervals)
    
    :param time: Time to run the simulation (in seconds)
    :type time: float
    :param N: Number of intervals of the metal bar
    :type N: int
    :param L: Length of the metal bar (in meters)
    :type L: float
    :param X0: Temperature at the beginning of the bar
    :type X0: int
    :param X1: Temperature at the end of the bar
    :type X1: int
    :param T0: Initial temperature of the bar (except the bar tips')
    :type T0: int

    :returns: matrix with rows = time, columns = position 
    :rtype: np.matrix
    '''

    dx = L/N

    global dt

    # rounds dt to get the most significant digit
    dt = round(0.2*dx*dx,-int(floor(log10(abs(0.2*dx*dx)))))

    number_time_interval = int(time/dt)

    temp_matrix = np.zeros((number_time_interval+1,N+1))

    # Boundary Conditions - temp_matrix[time,length]
    temp_matrix[0,:] = T0
    temp_matrix[:,0] = X0
    temp_matrix[:,N] = X1

    for k in range(0,number_time_interval):

        for i in range(1,N):
            
            temp_matrix[k+1,i] = temp_matrix[k,i] + (dt/(dx*dx)) *(temp_matrix[k,i+1] - 2*temp_matrix[k,i] + temp_matrix[k,i-1])

    return temp_matrix

def Q3(time,N,L,X0,X1)->np.matrix:
    ''' Returns the temperature of the bar on N intervals of the bar, 
    in different points of time (dt intervals)
    
    :param time: Time to run the simulation (in seconds)
    :type time: float
    :param N: Number of intervals of the metal bar
    :type N: int
    :param L: Length of the metal bar (in meters)
    :type L: float
    :param X0: Temperature at the beginning of the bar
    :type X0: int
    :param X1: Temperature at the end of the bar
    :type X1: int

    :returns: matrix with rows = time, columns = position 
    :rtype: np.matrix
    '''

    dx = L/N

    global dt

    # rounds dt to get the most significant digit
    dt = round(0.2*dx*dx,-int(floor(log10(abs(0.2*dx*dx)))))

    number_time_interval = int(time/dt)

    temp_matrix = np.zeros((number_time_interval+1,N+1))

    # Boundary Conditions - temp_matrix[time,length]
    temp_matrix[:,0] = X0
    temp_matrix[:,N] = X1

    for k in range(0,number_time_interval):

        for i in range(1,N):

            temp_matrix[0,i] = np.sin((np.pi/2)*(i*dx))
            
            temp_matrix[k+1,i] = temp_matrix[k,i] + (dt/(dx*dx))*(temp_matrix[k,i+1] - 2*temp_matrix[k,i] + temp_matrix[k,i-1])
        
    return temp_matrix

################################################################################
#                                                                              #
#                             MATHEMATICAL PROOFS                              #
#                                                                              #
################################################################################

def proof_Q1(time,N,L,X0,X1,T0,iterations)->np.matrix:
    ''' Returns the temperature of the bar on N intervals of the bar, 
    in different points of time (dt intervals)
    
    :param time: Time to run the simulation (in seconds)
    :type time: float
    :param N: Number of intervals of the metal bar
    :type N: int
    :param L: Length of the metal bar (in meters)
    :type L: float
    :param X0: Temperature at the beginning of the bar
    :type X0: int
    :param X1: Temperature at the end of the bar
    :type X1: int
    :param T0: Initial temperature of the bar (except the bar tips')
    :type T0: int
    :param iterations: number of iterations for the mathematical proof
    :type iterations: int

    :returns: matrix with rows = time, columns = position 
    :rtype: np.matrix
    '''

    dx = L/N

    # rounds dt to get the most significant digit
    dt = round(0.2*dx*dx,-int(floor(log10(abs(0.2*dx*dx)))))

    number_time_interval = int(time/dt)

    temp_matrix = np.zeros((number_time_interval+1,N+1))

    # Boundary Conditions - temp_matrix[time,length]
    temp_matrix[0,:] = T0
    temp_matrix[:,0] = X0
    temp_matrix[:,N] = X1

    iterations = iterations

    for k in range(0,number_time_interval):

        for i in range(1,N):

            for n in range(1,iterations):

                fraction = 4/((2*n-1)*np.pi)
                sin = np.sin((2*n-1)*np.pi*i*dx)
                exp = np.exp((-(2*n-1)**2)*(np.pi**2)*(k*dt))
            
                temp_matrix[k,i] += fraction*sin*exp

    return temp_matrix

def proof_Q2(time,N,L,X0,X1,T0,iterations)->np.matrix:
    ''' Returns the temperature of the bar on N intervals of the bar, 
    in different points of time (dt intervals)
    
    :param time: Time to run the simulation (in seconds)
    :type time: float
    :param N: Number of intervals of the metal bar
    :type N: int
    :param L: Length of the metal bar (in meters)
    :type L: float
    :param X0: Temperature at the beginning of the bar
    :type X0: int
    :param X1: Temperature at the end of the bar
    :type X1: int
    :param T0: Initial temperature of the bar (except the bar tips')
    :type T0: int
    :param iterations: number of iterations for the mathematical proof
    :type iterations: int

    :returns: matrix with rows = time, columns = position 
    :rtype: np.matrix
    '''

    dx = L/N

    # rounds dt to get the most significant digit
    dt = round(0.2*dx*dx,-int(floor(log10(abs(0.2*dx*dx)))))

    number_time_interval = int(time/dt)

    temp_matrix = np.zeros((number_time_interval+1,N+1))

    # Boundary Conditions - temp_matrix[time,length]
    temp_matrix[0,:] = T0
    temp_matrix[:,0] = X0
    temp_matrix[:,N] = X1

    iterations = iterations

    for k in range(0,number_time_interval):

        for i in range(1,N):

            sum = 0

            for n in range(1,iterations):

                fraction = 2/(n*np.pi)
                sin = np.sin((n)*np.pi*i*dx)
                exp = np.exp((-n**2)*(np.pi**2)*(k*dt))
            
                sum += fraction*sin*exp

                temp_matrix[k,i] = 1 - (i*dx) - sum

    return temp_matrix

def proof_Q3(time,N,L,X0,X1)->np.matrix:
    ''' Returns the temperature of the bar on N intervals of the bar, 
    in different points of time (dt intervals)
    
    :param time: Time to run the simulation (in seconds)
    :type time: float
    :param N: Number of intervals of the metal bar
    :type N: int
    :param L: Length of the metal bar (in meters)
    :type L: float
    :param X0: Temperature at the beginning of the bar
    :type X0: int
    :param X1: Temperature at the end of the bar
    :type X1: int

    :returns: matrix with rows = time, columns = position 
    :rtype: np.matrix
    '''

    dx = L/N

    # rounds dt to get the most significant digit
    dt = round(0.2*dx*dx,-int(floor(log10(abs(0.2*dx*dx)))))

    number_time_interval = int(time/dt)

    temp_matrix = np.zeros((number_time_interval+1,N+1))

    # Boundary Conditions - temp_matrix[time,length]
    temp_matrix[:,0] = X0
    temp_matrix[:,N] = X1

    for k in range(0,number_time_interval):

        for i in range(1,N):

            sin = np.sin((np.pi/2)*(i*dx))
            exp = np.exp(((-np.pi**2)*(k*dt))/4)

            temp_matrix[k,i] = exp*sin

    return temp_matrix

################################################################################
#                                                                              #
#                                  GRAPHS                                      #
#                                                                              #
################################################################################

def graph(time,proof,experimental,L,N,question):
    ''' Plots grah for proof and experimental data for the 3 questions
    
    :param time: Time to run the simulation (in seconds)
    :type time: float
    :param proof: Proof matrix
    :type proof: np.matrix
    :param experimental: Computed matrix
    :type experimental: np.matrix
    :param L: Length of the metal bar (in meters)
    :type L: float
    :param N: Number of intervals of the metal bar
    :type N: int
    :param question: Question name
    :type question: str
    
    :returns: 2x2 image
    :rtype: .png
    '''

    timeArray = np.arange(0,time,dt)

    # Time (in seconds) for each plotted phase
    # beginning, middle 1, middle 2, end
    timeList = [
        1,
        int(len(timeArray)/3),
        int((2*len(timeArray))/3),
        int(len(timeArray))-2
        ]

    positionArray = np.linspace(0,L,N+1)

    sns.set_theme(style="whitegrid")

    fig, axs = plt.subplots(figsize = (15,10),ncols = 2,nrows=2)

    dataColor = list(zip([experimental,proof],['#e63946','#219ebc']))

    for index,t in enumerate(timeList):

        col = index-2 if index >=2 else index

        for data,color in dataColor:

            sns.lineplot(
                x = positionArray,
                y = data[t,:],
                ax=axs[index//2,col],
                color = color,
                linewidth=2)

            axs[index//2,col].set_title(
                f'Time = {round(t*dt,2)} (seconds)',
                fontweight = 'bold',
                fontsize=17
                )

            axs[index//2,col].set_xlabel(
                'Position (m)',
                fontsize=15
                )

            axs[index//2,col].set_ylabel(
                'Temperature (°C)',
                fontsize=15
                )

    fig.suptitle(
        'Computed x Exact Solution',
        fontsize= 22,
        fontweight = 'bold',
        y = 1.01
    )    
    
    plt.legend(
        ['Computed', 'Exact Solution'],
        bbox_to_anchor=(-0.17, 2.60, 0.52, 0), 
        ncol=2,
        fontsize=18,
        frameon=False
        )
            
    # set the spacing between subplots
    plt.subplots_adjust( 
                    wspace=0.2, 
                    hspace=0.3
                )

    fig.savefig(f'./images/{question}.png',bbox_inches='tight')