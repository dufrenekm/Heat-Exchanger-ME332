# Heat Transfer Final Project
# Fall 2021
# Kyle and Marina

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

class heat_exchanger:
    def __init__(self):
        # Givens
        self.cp = 1e3  # J/(kg*K)
        self.t_c_in = 25  # C
        self.t_c_out = 225  # C
        self.t_h_in = 425  # C
        self.m_dot = 10  # kg/s
        self.h_given = 150  # W/(m^2*K)

        # Initial calculations for UA requirements
        self.c_c = self.m_dot*self.cp
        self.c_h = self.m_dot*self.cp
        self.c_min = min(self.c_c, self.c_h)
        self.c_max = max(self.c_c, self.c_h)
        self.c_r = self.c_min/self.c_max
        self.q = self.c_c*(self.t_c_out-self.t_c_in)
        self.q_max = self.c_min*(self.t_h_in-self.t_c_in)
        self.e = self.q/self.q_max
        self.ntu = -np.log(1+(1/self.c_r)*np.log(1-self.e*self.c_r))
        self.ua = self.c_min*self.ntu
        self.a_req = 2*self.ua/self.h_given

        # Properties cold
        self.k_c = 26.3e-3
        self.c_1 = .229
        self.c_2 = .97
        self.m = .632
        self.rho_c = 1.1614
        self.mu_c = 184.6e-7
        self.pr = .707

        # Properties hot, air at 700K
        self.k_h = 52.4e-3
        self.rho_h = .4975
        self.mu_h = 338.8e-7
        self.pr_h = .695

        # Initial spacing, we loop through all of these values
        self.diam = [.001, .01, .03, .05, .07, .1, .15, .2]
        self.num_tubes_x = [5, 10, 15, 20, 25, 30, 35, 100, 500, 1000]
        self.num_tubes_y = [5, 10, 15, 20, 25, 30, 35, 100, 500, 1000]
        self.lll = [.5, .7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 6.5, 7, 10, 15, 20]

    def update_spacing(self, i, j, k, q):
        # Calculate the spacing and dimensions of the heat exchanger
        st = 2*self.diam[i]
        sl = 2*self.diam[i]
        x = self.num_tubes_x[j]*st
        y = self.num_tubes_y[k]*sl
        return [st, sl, x, y]

    def cold_calculations(self, length, i, q):
        # Calculate h_bar for the cold stream
        v = self.m_dot/(self.rho_c*length[2]*self.lll[q])
        v_max = length[0]*v/(length[0]-self.diam[i])
        re_d = self.rho_c*v_max*self.diam[i]/self.mu_c
        nu_d = 1.13*self.c_1*self.c_2*re_d**self.m*self.pr**(1/3)
        h_bar = nu_d*self.k_c/self.diam[i]
        return h_bar

    def hot_calculations(self, i, j, k):
        # Calculate h_bar for the hot stream
        v_t = self.m_dot*4/(self.diam[i]**2*np.pi*self.rho_h*self.num_tubes_x[j]*self.num_tubes_y[k])
        re_d = self.rho_h*v_t*self.diam[i]/self.mu_h
        nu_d = .023*re_d**(4/5)*self.pr_h**.3
        h_bar = nu_d * self.k_h/self.diam[i]
        return h_bar

    def check_valid(self, i, j, k, q):
        # Check that the combination meets the requirement not dimension more than 2 times another
        retu = False
        spacin = self.update_spacing(i,j,k,q)
        if spacin[2] <= 2*spacin[3] and spacin[3] <= 2*spacin[2]:
            if self.lll[q] <= 2*spacin[3] and spacin[3] <= 2*self.lll[q]:
                if self.lll[q] <= 2*spacin[2] and spacin[2] <= 2*self.lll[q]:
                    retu = True
        return retu

    def new_area(self, min_h):
        return self.ua/min_h

    def inital_iteration(self, ):
        good_c = []
        good_h = []
        new_are = []
        min_h_array = []
        is_array = []
        js_array = []
        ks_array = []
        qs_array = []
        x_array = []
        y_array = []
        area_array = []

        # Iterate through the arrays of possible lengths, diamters, and numbers of tubes
        for i, diam in enumerate(self.diam):
            for j, num_x in enumerate(self.num_tubes_x):
                for k, num_y in enumerate(self.num_tubes_y):
                    for q, ll in enumerate(self.lll):
                        if self.check_valid(i, j, k, q):
                            lengths = self.update_spacing(i, j, k, q)
                            h_bar_cold = self.cold_calculations(lengths, i, q)
                            h_bar_hot = self.hot_calculations(i, j, k)

                            # Check that UA is close to the required value
                            new_ua = (1/(1/h_bar_cold+1/h_bar_hot))*np.pi*ll*diam*num_y*num_x
                            if  self.ua+1000> new_ua > (self.ua):
                                is_array.append(diam)
                                js_array.append(num_x)
                                ks_array.append(num_y)
                                qs_array.append(ll)
                                good_c.append(h_bar_cold)
                                good_h.append(h_bar_hot)
                                min_h_array.append(min(h_bar_hot, h_bar_cold))
                                new_are.append(new_ua)
                                x_array.append(lengths[0]*num_x)
                                y_array.append(lengths[1] * num_y)
                                area_array.append(np.pi*ll*diam*num_y*num_x)

        # Output a table of possible options to view
        fig = go.Figure(data=[go.Table(
            header=dict(values=['Area','new_ua', 'Min h_bar', 'H_c', 'H_h', 'Diameter', 'Num T', 'Num L', 'L', 'x', 'y'],
                        line_color='darkslategray',
                        fill_color='lightskyblue',
                        align='left'),
            cells=dict(values=[area_array,new_are, min_h_array, good_c, good_h, is_array, js_array, ks_array, qs_array, x_array, y_array],
                       line_color='darkslategray',
                       fill_color='lightcyan',
                       align='left'))
        ])

        fig.update_layout(width=1000, height=2000)
        fig.show()


if __name__ == "__main__":
    # This runs the code above to output the table of possible combinations
    heater = heat_exchanger()
    heater.inital_iteration()

    # This creates another instance of heat_exchanger() to be used to evaluate our choice with fouling
    second_heater = heat_exchanger()
    second_heater.diam = [.05]
    second_heater.num_tubes_x = [35]
    second_heater.num_tubes_y = [35]
    second_heater.lll = [6.75]
    lengths = second_heater.update_spacing(0, 0, 0, 0)
    h_bar_cold = second_heater.cold_calculations(lengths, 0, 0)
    h_bar_hot = second_heater.hot_calculations(0, 0, 0)
    print('H_h=', h_bar_hot, ', H_c=', h_bar_cold)
    fouling = .0004
    area = np.pi*second_heater.diam[0]*second_heater.lll[0]*second_heater.num_tubes_y[0]*second_heater.num_tubes_x[0]
    updated_u = 1/(2*.0004+1/(h_bar_hot)+1/(h_bar_cold))
    print('Updated Ua with fouling: ', updated_u*area)

    # Loop through T_c_o and generate area values for a scatterplot
    tco = np.linspace(225, 350, 100)
    q = second_heater.c_c * (tco - second_heater.t_c_in)
    e = q / second_heater.q_max
    ntu = np.zeros(100)
    for i in range(100):
        ntu[i] = -1*np.log(1+np.log(1-e[i]))
    ua = second_heater.c_min * ntu
    a_req = ua / updated_u
    # Generate the scatter plot and save it
    plt.scatter(tco, a_req)
    plt.xlabel("Temperature Cold Out (C)")
    plt.ylabel("Area Required (m^2)")
    plt.title("Heat Exchanger Effectiveness")
    plt.savefig('ht322_final_project.png')
    plt.show()







