import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

class Lokta_Volterra_normal():
    def __init__(self,x0, y0, r=1, a=1, b=1, d=1, s=0):
        self.x0 = x0
        self.y0 = y0
        self.r = r
        self.a = a
        self.b = b
        self.d = d
        self.s = s


    def solutions(self,T, step):
 
        time = np.arange(0, T, step)

        x = np.zeros(len(time))
        y = np.zeros(len(time))
    
        x[0] = self.x0
        y[0] = self.y0
        for i in range(1,len(time)):
            x[i] = x[i-1] + (x[i-1]*(self.r - self.a*y[i-1] - x[i-1]*self.s))*step
            y[i] = y[i-1] + (y[i-1]*(self.b*x[i-1] - self.d) )*step 
        return x,y,time
        
    def plot_solutions(self,title = '',T=50, step=0.0001):
        
        x_sol,y_sol,time = self.solutions(T=T, step=step)
                             
        fig, ax = plt.subplots(2,figsize=(7,6))
        fig.suptitle(title)
        ax[0].plot(time, x_sol, label = 'Prey')
        ax[0].plot(time, y_sol, label = 'Predator')
        ax[0].legend()
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel('Population Concentration')
    
        ax[1].plot(x_sol, y_sol,alpha=0.5,c='black')
        ax[1].arrow(x_sol[0], y_sol[0],x_sol[20]-x_sol[0],y_sol[20]-y_sol[0] ,shape='full', color='black', lw=0, head_width=.05)
        ax[1].annotate('x0, y0', xy=(x_sol[0], y_sol[0]), xytext=(-40, -10), textcoords='offset points', c='black')
        
        stable_point = ( self.d/self.b, self.r/self.a - self.d*self.s/(self.a*self.b) )
        ax[1].scatter(stable_point[0], stable_point[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x}$, $\hat{y}$', xy=stable_point, xytext=(-8, 10), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        plt.tight_layout()

class Lokta_Volterra_two_prey():
    def __init__(self,x1_0, x2_0, y0, r=1, a=1, b=1, d=1, eps=(0,0), rho=1/2, change=1, type = 'change'):
        self.x1_0 = x1_0
        self.x2_0 = x2_0
        self.y0 = y0
        self.r = r
        self.a = a
        self.b = b
        self.d = d
        self.eps = eps
        self.rho = rho
        self.change = change
        self.type = type

    def solutions(self,T, step):
 
        time = np.arange(0, T, step)

        x1 = np.zeros(len(time))
        x2 = np.zeros(len(time))
        y = np.zeros(len(time))
    
        x1[0] = self.x1_0
        x2[0] = self.x2_0
        y[0] = self.y0

        if self.type=='change':
            eps_list = [self.eps[1]]*int(len(time))
    
            eps_list[0:int(len(time)*self.change)]= [self.eps[0]]*int(len(time)*self.change)
        else:
            eps_list = self.eps
            
        for i in range(1,len(time)):
            x1[i] = x1[i-1] + (self.r*(self.rho*x1[i-1] + (1-self.rho)*x2[i-1]) - self.a*(1+eps_list[i])*x1[i-1]*y[i-1])*step
            x2[i] = x2[i-1] + (self.r*((1-self.rho)*x1[i-1] + self.rho*x2[i-1]) - self.a*(1-eps_list[i])*x2[i-1]*y[i-1])*step
            y[i] = y[i-1] + (y[i-1]*(self.b*((1+eps_list[i])*x1[i-1] + (1-eps_list[i])*x2[i-1]) - self.d) )*step
        return x1, x2, y, time
        
    def equa(self,p):
        x1, x2, y = p
        r = self.r
        a = self.a
        b = self.b
        d = self.d
        rho = self.rho
        eps = self.eps[0]
        
        x1_dot = r*(rho*x1+(1-rho)*x2 - a*(1+eps)*x1*y)
        x2_dot = r*(rho*x2+(1-rho)*x1 - a*(1-eps)*x2*y)
        y_dot = y*(b*((1+eps)*x1+(1-eps)*x2)-d)
        return (x1_dot,x2_dot,y_dot)

        
    def plot_solutions(self,title = '',T=50, step=0.001, line_change=False):
        
        x1_sol, x2_sol, y_sol, time = self.solutions(T=T, step=step)

        solution_stable = fsolve(self.equa, (1, 1, 1))
        print(solution_stable)
                             
        fig, ax = plt.subplots(2,figsize=(7,6))
        fig.suptitle(title)
        ax[0].plot(time, x1_sol, label = 'Prey 1')
        ax[0].plot(time, x2_sol, label = 'Prey 2')
        ax[0].plot(time, y_sol, label = 'Predator')
        ax[0].legend()
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel('Population Concentration')

        # PREY 1
    
        ax[1].plot(x1_sol, y_sol,alpha=0.5, label = 'Prey 1')
        ax[1].arrow(x1_sol[0], y_sol[0],x1_sol[10]-x1_sol[0],y_sol[10]-y_sol[0] ,shape='full', color='black',length_includes_head=True, lw=0, head_width=.05)
        ax[1].annotate(r'$x1_0, y_0$', xy=(x1_sol[0], y_sol[0]), xytext=(0, +10), textcoords='offset points', c='black')

        stable_point_1 = ( solution_stable[0] , solution_stable[2] )
        ax[1].scatter(stable_point_1[0], stable_point_1[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x_1}$, $\hat{y}$', xy=stable_point_1, xytext=(-8, 10), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        plt.tight_layout()

        # PREY 2
                          
        ax[1].plot(x2_sol, y_sol,alpha=0.5,c='orange',label='Prey 2')
        ax[1].arrow(x2_sol[0], y_sol[0],x2_sol[10]-x2_sol[0],y_sol[10]-y_sol[0] ,shape='full', color='black',length_includes_head=True, lw=0, head_width=.05)
        ax[1].annotate(r'$x2_0, y_0$', xy=(x2_sol[0], y_sol[0]), xytext=(-40, -10), textcoords='offset points', c='black')
        
        stable_point_2 = ( solution_stable[1] , solution_stable[2] )
        ax[1].scatter(stable_point_2[0], stable_point_2[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x_2}$, $\hat{y}$', xy=stable_point_2, xytext=(-8, 10), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        ax[1].legend()
        if line_change == True:
            ax[0].axvline(int(len(time)*self.change*step),ls=':', alpha = 0.5, c='red')
        
        plt.tight_layout()
        
    def max_amplitude(self,position=100,step = 0.001):

        x1_sol, x2_sol, y_sol, time = self.solutions(T=position+20, step=step)

        lim_inf = int((position-10)/step)
        lim_sup = int((position+10)/step)
        x1_max = np.max(x1_sol[lim_inf:lim_sup])
        x2_max = np.max(x2_sol[lim_inf:lim_sup])
        y_max = np.max(y_sol[lim_inf:lim_sup])
        return x1_max, x2_max, y_max

    def size_perturbations(self, initial, each, T=1000, step = 0.001):
        
        x1_sol, x2_sol, y_sol, time = self.solutions(T=T, step=step)

        for i in range(int((T-initial)/each)):
            x1_max = np.max(x1_sol[initial+i*int(each/step):initial+(i+1)*int(each/step)])
            x2_max = np.max(x2_sol[initial+i*int(each/step):initial+(i+1)*int(each/step)])
            y_max = np.max(y_sol[initial+i*int(each/step):initial+(i+1)*int(each/step)])

            x1_max = x1_max * 1/int((T-initial)/each)
            x2_max = x2_max * 1/int((T-initial)/each)
            y_max = y_max * 1/int((T-initial)/each)
        return(x1_max,x2_max,y_max)
            
            
        
class Lokta_Volterra_two_predator():
    def __init__(self,x_0, y1_0, y2_0, r=1, a=1, b=1, d=1, delta=(0,0), rho=1/2, change=1/2, type='change'):
        self.x_0 = x_0
        self.y1_0 = y1_0
        self.y2_0 = y2_0
        self.r = r
        self.a = a
        self.b = b
        self.d = d
        self.delta = delta
        self.rho = rho
        self.change = change
        self.type = type
    def solutions(self,T, step):

        time = np.arange(0, T, step)

        x = np.zeros(len(time))
        y1 = np.zeros(len(time))
        y2 = np.zeros(len(time))

        if self.type == 'change':
            delta_list = [self.delta[1]]*int(len(time))
    
            delta_list[0:int(len(time)*self.change)]= [self.delta[0]]*int(len(time)*self.change)
        else:
            delta_list = self.delta
    
        x[0] = self.x_0
        y1[0] = self.y1_0
        y2[0] = self.y2_0
        for i in range(1,len(time)):
            x[i] = x[i-1] + x[i-1]*(self.r - self.a*(y1[i-1]+y2[i-1]))*step
            y1[i] = y1[i-1] + (-self.d*(1+delta_list[i])*y1[i-1] + self.b*(self.rho*y1[i-1]+(1-self.rho)*y2[i-1])*x[i-1])*step 
            y2[i] = y2[i-1] + (-self.d*(1-delta_list[i])*y2[i-1] + self.b*((1-self.rho)*y1[i-1]+self.rho*y2[i-1])*x[i-1])*step 
        return x, y1, y2, time
        
    def equa(self,p):
        x, y1, y2 = p
        r = self.r
        a = self.a
        b = self.b
        d = self.d
        rho = self.rho
        delta = self.delta[0]
        x_dot = x*(r-a*(y1+y2))
        y1_dot = -d*(1+delta)*y1 + b*(rho*y1 + (1-rho)*y2)*x
        y2_dot = -d*(1-delta)*y2 + b*(rho*y2 + (1-rho)*y1)*x
        return (x_dot,y1_dot,y2_dot)

    def plot_solutions(self,title = '',T=50, step=0.001, line_change=False):
        
        x_sol, y1_sol, y2_sol, time = self.solutions(T=T, step=step)

        solution_stable = fsolve(self.equa, (1, 1, 1))
        print(solution_stable)
                             
        fig, ax = plt.subplots(2,figsize=(7,6))
        fig.suptitle(title)
        ax[0].plot(time, x_sol, label = 'Prey')
        ax[0].plot(time, y1_sol, label = 'Predator 1')
        ax[0].plot(time, y2_sol, label = 'Predator 2')
        ax[0].legend()
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel('Population Concentration')

        # PREY 1
    
        ax[1].plot(x_sol, y1_sol,alpha=0.5,c='orange', label = 'Predator 1')
        ax[1].arrow(x_sol[0], y1_sol[0],x_sol[10]-x_sol[0],y1_sol[10]-y1_sol[0] ,shape='full', color='black',length_includes_head=True, lw=0, head_width=.05)
        ax[1].annotate(r'$x_0, y1_0$', xy=(x_sol[0], y1_sol[0]), xytext=(0, +10), textcoords='offset points', c='black')

        
        stable_point_1 = ( solution_stable[0] , solution_stable[1] )
        ax[1].scatter(stable_point_1[0], stable_point_1[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x}$, $\hat{y_1}$', xy=stable_point_1, xytext=(+10, -10), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        plt.tight_layout()

        # PREY 2
                          
        ax[1].plot(x_sol, y2_sol,alpha=0.5,c='green',label='Predator 2')
        ax[1].arrow(x_sol[0], y2_sol[0],x_sol[10]-x_sol[0],y2_sol[10]-y2_sol[0] ,shape='full', color='black',length_includes_head=True, lw=0, head_width=.05)
        ax[1].annotate(r'$x_0, y2_0$', xy=(x_sol[0], y2_sol[0]), xytext=(-40, -10), textcoords='offset points', c='black')
        
        stable_point_2 = ( solution_stable[0] , solution_stable[2] )
        ax[1].scatter(stable_point_2[0], stable_point_2[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x}$, $\hat{y_2}$', xy=stable_point_2, xytext=(+10, 0), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        ax[1].legend()

        if line_change == True:
            ax[0].axvline(int(len(time)*self.change*step),ls=':', alpha = 0.5, c='red')
        
        plt.tight_layout()

    def max_amplitude(self,position=100,step = 0.001):

        x_sol, y1_sol, y2_sol, time = self.solutions(T=position+20, step=step)

        lim_inf = int((position-10)/step)
        lim_sup = int((position+10)/step)
        x_max = np.max(x_sol[lim_inf:lim_sup])
        y1_max = np.max(y1_sol[lim_inf:lim_sup])
        y2_max = np.max(y2_sol[lim_inf:lim_sup])
        return x_max, y1_max, y2_max

    def size_perturbations(self, initial, each, T=1000, step = 0.001):
        
        x1_sol, x2_sol, y_sol, time = self.solutions(T=T, step=step)

        for i in range(int((T-initial)/each)):
            x1_max = np.max(x1_sol[initial+int(i*each)/step:initial+(i+1)*int(each/step)])
            x2_max = np.max(x2_sol[initial+int(i*each)/step:initial+(i+1)*int(each/step)])
            y_max = np.max(y_sol[initial+int(i*each)/step:initial+(i+1)*int(each/step)])

            x1_max = x1_max * 1/int((T-initial)/each)
            x2_max = x2_max * 1/int((T-initial)/each)
            y_max = y_max * 1/int((T-initial)/each)
        return(x1_max,x2_max,y_max)
            
class Lokta_Volterra_normal_simulation():
    def __init__(self,x0, y0, r=1, a=1, b=1, d=1, s=0,N=100):
        self.x0 = x0
        self.y0 = y0
        self.r = r
        self.a = a
        self.b = b
        self.d = d
        self.s = s
        self.N = N


    def solutions(self,T, step):
        r = self.r
        a = self.a
        b = self.b
        d = self.d
        s = self.s
        N = self.N
        time = np.arange(0, T, step)

        x = np.zeros(len(time))
        y = np.zeros(len(time))
    
        x[0] = self.x0
        y[0] = self.y0
        for i in range(1,len(time)):
            #reproduction prey
            rand_r = np.random.random(int(x[i-1]*N))
            new_r = np.sum((rand_r<(r*step)).astype(int))/N
            x[i] = x[i-1] + new_r
            
            #death prey
            rand_a = np.random.random(int(x[i-1]*N))
            new_a = np.sum((rand_a<(a*step*y[i-1])).astype(int))/N
            x[i] = x[i] - new_a
        
            #carrying capacity prey
            rand_s = np.random.random(int(x[i-1]*N))
            new_s= np.sum((rand_s<(s*step*x[i-1])).astype(int))/N
            x[i] = x[i] - new_s
        
            #reproduction predator
            rand_b = np.random.random(int(y[i-1]*N))
            new_b = np.sum((rand_b<(b*step*x[i-1])).astype(int))/N
            y[i] = y[i-1] + new_b
            
        
            #death predator
            rand_d = np.random.random(int(y[i-1]*N))
            new_d = np.sum((rand_d<(d*step)).astype(int))/N
            y[i] = y[i] - new_d
        return x, y, time
        
    def plot_solutions(self,title = '',T=50, step=0.001):
        
        x_sol,y_sol,time = self.solutions(T=T, step=step)
                             
        fig, ax = plt.subplots(2,figsize=(7,6))
        fig.suptitle(title)
        ax[0].plot(time, x_sol, label = 'Prey')
        ax[0].plot(time, y_sol, label = 'Predator')
        ax[0].legend()
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel('Population Concentration')
    
        ax[1].plot(x_sol, y_sol,alpha=0.5,c='black')
        ax[1].arrow(x_sol[0], y_sol[0],x_sol[20]-x_sol[0],y_sol[20]-y_sol[0] ,shape='full', color='black', lw=0, head_width=.05)
        ax[1].annotate('x0, y0', xy=(x_sol[0], y_sol[0]), xytext=(-40, -10), textcoords='offset points', c='black')
        
        stable_point = ( self.d/self.b, self.r/self.a - self.d*self.s/(self.a*self.b) )
        ax[1].scatter(stable_point[0], stable_point[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x}$, $\hat{y}$', xy=stable_point, xytext=(-8, 10), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        plt.tight_layout()

class Lokta_Volterra_two_prey_simulation():
    def __init__(self,x1_0, x2_0, y0, r=1, a=1, b=1, d=1, eps=(0,0), rho=1/2, N=100, change=1,type='change'):
        self.x1_0 = x1_0
        self.x2_0 = x2_0
        self.y0 = y0
        self.r = r
        self.a = a
        self.b = b
        self.d = d
        self.eps = eps
        self.N = N
        self.rho = rho
        self.change = change
        self.type = type


    def solutions(self,T, step):
        
        r = self.r
        a = self.a
        b = self.b
        d = self.d
        N = self.N
        eps = self.eps
        time = np.arange(0, T, step)

        x1 = np.zeros(len(time))
        x2 = np.zeros(len(time))
        y = np.zeros(len(time))
    
        x1[0] = self.x1_0
        x2[0] = self.x2_0
        y[0] = self.y0
        
        if self.type == 'change':
            eps_list = [self.eps[1]]*int(len(time))
    
            eps_list[0:int(len(time)*self.change)]= [self.eps[0]]*int(len(time)*self.change)
        else:
            eps_list = self.eps
        
        for i in range(1,len(time)):
            #reproduction prey 1 to 1
            rand_r11 = np.random.random(int(x1[i-1]*N))
            new_r11 = np.sum((rand_r11<(r*self.rho*step)).astype(int))/N
            x1[i] = x1[i-1] + new_r11
            
            
            #reproduction prey 1 to 2
            rand_r12 = np.random.random(int(x1[i-1]*N))
            new_r12 = np.sum((rand_r12<(r*(1-self.rho)*step)).astype(int))/N
            x2[i] = x2[i-1] + new_r12

            #reproduction prey 2 to 1
            rand_r21 = np.random.random(int(x2[i-1]*N))
            new_r21 = np.sum((rand_r21<(r*(1-self.rho)*step)).astype(int))/N
            x1[i] = x1[i] + new_r21
            
            
            #reproduction prey 2 to 2
            rand_r22 = np.random.random(int(x2[i-1]*N))
            new_r22 = np.sum((rand_r22<(r*self.rho*step)).astype(int))/N
            x2[i] = x2[i] + new_r22

            
            #death prey 1
            rand_a1 = np.random.random(int(x1[i-1]*N))
            new_a1 = np.sum((rand_a1<(a*(1+eps_list[i])*step*y[i-1])).astype(int))/N
            x1[i] = x1[i] - new_a1

            #death prey 2
            rand_a2 = np.random.random(int(x2[i-1]*N))
            new_a2 = np.sum((rand_a2<(a*(1-eps_list[i])*step*y[i-1])).astype(int))/N
            x2[i] = x2[i] - new_a2
        
            #reproduction predator
            rand_b = np.random.random(int(y[i-1]*N))
            new_b = np.sum((rand_b<(b*step*((1+eps_list[i])*x1[i-1]+(1-eps_list[i])*x2[i-1]))).astype(int))/N
            y[i] = y[i-1] + new_b
            
        
            #death predator
            rand_d = np.random.random(int(y[i-1]*N))
            new_d = np.sum((rand_d<(d*step)).astype(int))/N
            y[i] = y[i] - new_d
            
        return x1, x2, y, time
        
    def equa(self,p):
        x1, x2, y = p
        r = self.r
        a = self.a
        b = self.b
        d = self.d
        rho = self.rho
        eps = self.eps[0]
        
        x1_dot = r*(rho*x1+(1-rho)*x2 - a*(1+eps)*x1*y)
        x2_dot = r*(rho*x2+(1-rho)*x1 - a*(1-eps)*x2*y)
        y_dot = y*(b*((1+eps)*x1+(1-eps)*x2)-d)
        return (x1_dot,x2_dot,y_dot)

        
    def plot_solutions(self,title = '',T=50, step=0.001, line_change=False):
        
        x1_sol, x2_sol, y_sol, time = self.solutions(T=T, step=step)

        solution_stable = fsolve(self.equa, (1, 1, 1))
        print(solution_stable)
                             
        fig, ax = plt.subplots(2,figsize=(7,6))
        fig.suptitle(title)
        ax[0].plot(time, x1_sol, label = 'Prey 1')
        ax[0].plot(time, x2_sol, label = 'Prey 2')
        ax[0].plot(time, y_sol, label = 'Predator')
        ax[0].legend()
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel('Population Concentration')

        # PREY 1
    
        ax[1].plot(x1_sol, y_sol,alpha=0.5, label = 'Prey 1')
        ax[1].arrow(x1_sol[0], y_sol[0],x1_sol[2]-x1_sol[0],y_sol[2]-y_sol[0] ,shape='full', color='black',length_includes_head=True, lw=0, head_width=.05)
        ax[1].annotate(r'$x1_0, y_0$', xy=(x1_sol[0], y_sol[0]), xytext=(0, +10), textcoords='offset points', c='black')

        stable_point_1 = ( solution_stable[0] , solution_stable[2] )
        ax[1].scatter(stable_point_1[0], stable_point_1[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x_1}$, $\hat{y}$', xy=stable_point_1, xytext=(-8, 10), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        plt.tight_layout()

        # PREY 2
                          
        ax[1].plot(x2_sol, y_sol,alpha=0.5,c='orange',label='Prey 2')
        ax[1].arrow(x2_sol[0], y_sol[0],x2_sol[2]-x2_sol[0],y_sol[2]-y_sol[0] ,shape='full', color='black',length_includes_head=True, lw=0, head_width=.05)
        ax[1].annotate(r'$x2_0, y_0$', xy=(x2_sol[0], y_sol[0]), xytext=(-40, -10), textcoords='offset points', c='black')
        
        stable_point_2 = ( solution_stable[1] , solution_stable[2] )
        ax[1].scatter(stable_point_2[0], stable_point_2[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x_2}$, $\hat{y}$', xy=stable_point_2, xytext=(-8, 10), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        ax[1].legend()
        if line_change == True:
            ax[0].axvline(int(len(time)*self.change*step),ls=':', alpha = 0.5, c='red')
        
        plt.tight_layout()


    def size_perturbations(self, initial, each, T=1000, step = 0.001):
        
        x1_sol, x2_sol, y_sol, time = self.solutions(T=T, step=step)

        for i in range(int((T-initial)/each)):
            x1_max = np.max(x1_sol[initial+int(i*each)/step:initial+(i+1)*int(each/step)])
            x2_max = np.max(x2_sol[initial+int(i*each)/step:initial+(i+1)*int(each/step)])
            y_max = np.max(y_sol[initial+int(i*each)/step:initial+(i+1)*int(each/step)])

            x1_max = x1_max * 1/int((T-initial)/each)
            x2_max = x2_max * 1/int((T-initial)/each)
            y_max = y_max * 1/int((T-initial)/each)
        return(x1_max,x2_max,y_max)
        
    
        
class Lokta_Volterra_two_predator_simulation():
    def __init__(self,x_0, y1_0, y2_0, r=1, a=1, b=1, d=1, delta=(0,0), rho=1/2, change = 1/2, N=100, type = 'change'):
        self.x_0 = x_0
        self.y1_0 = y1_0
        self.y2_0 = y2_0
        self.r = r
        self.a = a
        self.b = b
        self.d = d
        self.delta = delta
        self.N = N
        self.rho = rho
        self.change = change
        self.type = type

    def solutions(self,T, step):

        r = self.r
        a = self.a
        b = self.b
        d = self.d
        N = self.N
        rho = self.rho

        time = np.arange(0, T, step)
        
        if self.type == 'change':
            delta_list = [self.delta[1]]*int(len(time))
    
            delta_list[0:int(len(time)*self.change)]= [self.delta[0]]*int(len(time)*self.change)
        else:
            delta_list = self.delta
        x = np.zeros(len(time))
        y1 = np.zeros(len(time))
        y2 = np.zeros(len(time))
    
        x[0] = self.x_0
        y1[0] = self.y1_0
        y2[0] = self.y2_0
        for i in range(1,len(time)):
            #reproduction prey
            rand_r = np.random.random(int(x[i-1]*N))
            new_r = np.sum((rand_r<(r*step)).astype(int))/N
            x[i] = x[i-1] + new_r
            
            #death prey 1
            rand_a = np.random.random(int(x[i-1]*N))
            new_a = np.sum((rand_a<(a*step*(y1[i-1]+y2[i-1]))).astype(int))/N
            x[i] = x[i] - new_a


            #reproduction predator 1 to 1
            rand_b11 = np.random.random(int(y1[i-1]*N))
            new_b11 = np.sum((rand_b11<(b*rho*step*x[i-1])).astype(int))/N
            y1[i] = y1[i-1] + new_b11
            
            #reproduction predator 1 to 2
            rand_b12 = np.random.random(int(y1[i-1]*N))
            new_b12 = np.sum((rand_b12<(b*(1-rho)*step*x[i-1])).astype(int))/N
            y2[i] = y2[i-1] + new_b12
            
            #reproduction predator 2 to 2
            rand_b2 = np.random.random(int(y2[i-1]*N))
            new_b2 = np.sum((rand_b2<(b*rho*step*x[i-1])).astype(int))/N
            y2[i] = y2[i] + new_b2
            
            #reproduction predator 2 to 1
            rand_b21 = np.random.random(int(y2[i-1]*N))
            new_b21 = np.sum((rand_b21<(b*(1-rho)*step*x[i-1])).astype(int))/N
            y1[i] = y1[i] + new_b21
            
        
            #death predator 1
            rand_d1 = np.random.random(int(y1[i-1]*N))
            new_d1 = np.sum((rand_d1<(d*(1+delta_list[i])*step)).astype(int))/N
            y1[i] = y1[i] - new_d1
            
            #death predator 2
            rand_d2 = np.random.random(int(y2[i-1]*N))
            new_d2 = np.sum((rand_d2<(d*(1-delta_list[i])*step)).astype(int))/N
            y2[i] = y2[i] - new_d2
            
        return x, y1, y2, time

    def equa(self,p):
        x, y1, y2 = p
        r = self.r
        a = self.a
        b = self.b
        d = self.d
        rho = self.rho
        delta = self.delta[0]
        x_dot = x*(r-a*(y1+y2))
        y1_dot = -d*(1+delta)*y1 + b*(rho*y1 + (1-rho)*y2)*x
        y2_dot = -d*(1-delta)*y2 + b*(rho*y2 + (1-rho)*y1)*x
        return (x_dot,y1_dot,y2_dot)
        
    def plot_solutions(self,title = '',T=50, step=0.001, line_change=False):
        
        x_sol, y1_sol, y2_sol, time = self.solutions(T=T, step=step)

        solution_stable = fsolve(self.equa, (1, 1, 1))
        print(solution_stable)
                             
        fig, ax = plt.subplots(2,figsize=(7,6))
        fig.suptitle(title)
        ax[0].plot(time, x_sol, label = 'Prey')
        ax[0].plot(time, y1_sol, label = 'Predator 1')
        ax[0].plot(time, y2_sol, label = 'Predator 2')
        ax[0].legend()
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel('Population Concentration')

        # PREDATOR 1
    
        ax[1].plot(x_sol, y1_sol,alpha=0.5,c='orange', label = 'Predator 1')
        ax[1].arrow(x_sol[0], y1_sol[0],x_sol[2]-x_sol[0],y1_sol[2]-y1_sol[0] ,shape='full', color='black',length_includes_head=True, lw=0, head_width=.05)
        ax[1].annotate(r'$x_0, y1_0$', xy=(x_sol[0], y1_sol[0]), xytext=(0, +10), textcoords='offset points', c='black')

        stable_point_1 = ( solution_stable[0] , solution_stable[1] )
        ax[1].scatter(stable_point_1[0], stable_point_1[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x}$, $\hat{y_1}$', xy=stable_point_1, xytext=(+10, -10), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        plt.tight_layout()

        # PREDATOR 2
                          
        ax[1].plot(x_sol, y2_sol,alpha=0.5,c='green',label='Predator 2')
        ax[1].arrow(x_sol[0], y2_sol[0],x_sol[2]-x_sol[0],y2_sol[2]-y2_sol[0] ,shape='full', color='black',length_includes_head=True, lw=0, head_width=.05)
        ax[1].annotate(r'$x_0, y2_0$', xy=(x_sol[0], y2_sol[0]), xytext=(-40, -10), textcoords='offset points', c='black')
        
        stable_point_2 = ( solution_stable[0] , solution_stable[2] )
        ax[1].scatter(stable_point_2[0], stable_point_2[1], marker='*', c='red')
        ax[1].annotate(r'$\hat{x}$, $\hat{y_2}$', xy=stable_point_2, xytext=(+10, 0), textcoords='offset points',c='red')
        ax[1].set_xlabel('Concentration prey')
        ax[1].set_ylabel('Concentration predator')
        ax[1].legend()
        if line_change == True:
            ax[0].axvline(int(len(time)*self.change*step),ls=':', alpha = 0.5, c='red')
        
        plt.tight_layout()


    def size_perturbations(self, initial, each, T=1000, step = 0.001):
        
        x1_sol, x2_sol, y_sol, time = self.solutions(T=T, step=step)

        for i in range(int((T-initial)/each)):
            x1_max = np.max(x1_sol[initial+int(i*each)/step:initial+(i+1)*int(each/step)])
            x2_max = np.max(x2_sol[initial+int(i*each)/step:initial+(i+1)*int(each/step)])
            y_max = np.max(y_sol[initial+int(i*each)/step:initial+(i+1)*int(each/step)])

            x1_max = x1_max * 1/int((T-initial)/each)
            x2_max = x2_max * 1/int((T-initial)/each)
            y_max = y_max * 1/int((T-initial)/each)
        return(x1_max,x2_max,y_max)

def maximum_plot_prey(pos=100,rho_step=0.01):
    rho_list = np.arange(0,1+rho_step,rho_step)
    x1_max = np.zeros(len(rho_list))
    x2_max = np.zeros(len(rho_list))
    y_max = np.zeros(len(rho_list))
    for i in range(len(rho_list)):
        x1_max[i] , x2_max[i], y_max[i] = Lokta_Volterra_two_prey(0.4,0.4,0.5,eps=(1/4,1/4),rho=rho_list[i]).max_amplitude(position=pos)
    plt.plot(rho_list, x1_max, label = 'Prey 1')
    plt.plot(rho_list, x2_max, label = 'Prey 2')
    plt.plot(rho_list, y_max, label = 'Predator')
    plt.legend()
    plt.xlabel(r'$\rho$')
    plt.ylabel(f'max value around {pos} ')
    plt.show()

def maximum_plot_predator(pos=100,rho_step=0.01):
    rho_list = np.arange(0,1+rho_step,rho_step)
    x_max = np.zeros(len(rho_list))
    y1_max = np.zeros(len(rho_list))
    y2_max = np.zeros(len(rho_list))
    for i in range(len(rho_list)):
        x_max[i] , y1_max[i], y2_max[i] = Lokta_Volterra_two_predator(0.7,0.2,0.2,delta=(1/3,+1/3),rho=rho_list[i],change = 1).max_amplitude(position=pos)
    plt.plot(rho_list, x_max, label = 'Prey')
    plt.plot(rho_list, y1_max, label = 'Predator 1')
    plt.plot(rho_list, y2_max, label = 'Predator 2')
    plt.legend()
    plt.xlabel(r'$\rho$')
    plt.ylabel(f'max value around {pos} ')
    plt.show()

def amplitude_plot_prey(rho_step=0.01,initial = 100, each = 10):
    var = 1/20
    step_ = 0.001
    epsilon_list = np.random.normal(1/4, var,int(1100))
    epsilon_list = np.repeat(epsilon_list, int(1/step_))
    rho_list = np.arange(0,0.95,rho_step)
    x1_max = np.zeros(len(rho_list))
    x2_max = np.zeros(len(rho_list))
    y_max = np.zeros(len(rho_list))
    for i in range(len(rho_list)):
        x1_max[i] , x2_max[i], y_max[i] = Lokta_Volterra_two_prey(2,2,3,eps=epsilon_list,rho=rho_list[i],change = 1, type='').size_perturbations(initial=initial,each=each)
    plt.plot(rho_list, x1_max, label = 'Prey 1')
    plt.plot(rho_list, x2_max, label = 'Prey 2')
    plt.plot(rho_list, y_max, label = 'Predator')
    plt.legend()
    plt.xlabel(r'$\rho$')
    plt.ylabel(f'aproximate of amplitude of oscillations ')
    plt.show()
    
def amplitude_plot_predator(rho_step=0.01,initial = 100, each = 10):
    var = 1/20
    step_ = 0.001
    epsilon_list = np.random.normal(1/4, var,int(1100))
    epsilon_list = np.repeat(epsilon_list, int(1/step_))
    rho_list = np.arange(0,0.95,rho_step)
    x1_max = np.zeros(len(rho_list))
    x2_max = np.zeros(len(rho_list))
    y_max = np.zeros(len(rho_list))
    for i in range(len(rho_list)):
        x1_max[i] , x2_max[i], y_max[i] = Lokta_Volterra_two_predator(2,3,3,delta=epsilon_list,rho=rho_list[i],change = 1, type='').size_perturbations(initial=initial,each=each)
    plt.plot(rho_list, x1_max, label = 'Prey')
    plt.plot(rho_list, x2_max, label = 'Predator 1')
    plt.plot(rho_list, y_max, label = 'Predator 2')
    plt.legend()
    plt.xlabel(r'$\rho$')
    plt.ylabel(f'aproximate of amplitude of oscillations ')
    plt.show()

def amplitude_plot_prey_simulation(rho_step=0.01,initial = 100, each = 10):
    step_ = 0.05
    rho_list = np.arange(0,0.95,rho_step)
    x1_max = np.zeros(len(rho_list))
    x2_max = np.zeros(len(rho_list))
    y_max = np.zeros(len(rho_list))
    for i in range(len(rho_list)):
        x1_max[i] , x2_max[i], y_max[i] = Lokta_Volterra_two_prey_simulation(2,2,3,eps=(1/4,1/4),rho=rho_list[i],change = 1, N=1000, type='').size_perturbations(initial=initial,each=each,T=500,step=step_)
    plt.plot(rho_list, x1_max, label = 'Prey 1')
    plt.plot(rho_list, x2_max, label = 'Prey 2')
    plt.plot(rho_list, y_max, label = 'Predator')
    plt.legend()
    plt.xlabel(r'$\rho$')
    plt.ylabel(f'aproximate of amplitude of oscillations ')
    plt.show()
    
def amplitude_plot_predator_simulation(rho_step=0.01,initial = 100, each = 10):
    var = 1/20
    step_ = 0.05
    rho_list = np.arange(0,0.95,rho_step)
    x1_max = np.zeros(len(rho_list))
    x2_max = np.zeros(len(rho_list))
    y_max = np.zeros(len(rho_list))
    for i in range(len(rho_list)):
        x1_max[i] , x2_max[i], y_max[i] = Lokta_Volterra_two_predator_simulation(2,3,3,delta=(1/3,1/3),rho=rho_list[i],change = 1, N=1000, type='').size_perturbations(initial=initial,each=each)
    plt.plot(rho_list, x1_max, label = 'Prey')
    plt.plot(rho_list, x2_max, label = 'Predator 1')
    plt.plot(rho_list, y_max, label = 'Predator 2')
    plt.legend()
    plt.xlabel(r'$\rho$')
    plt.ylabel(f'aproximate of amplitude of oscillations ')
    plt.show()