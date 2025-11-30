import streamlit as st
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sympy import symbols, simplify, Eq, diff, I, sin, plotting, sqrt, pi, exp, lambdify

tempo = st.selectbox('Selecione se é dependente do tempo ou se é independente: ', ['Dependente','Independente'])

dimensão = st.selectbox('Selecione o número de seleções', ['1D','2D','3D'])

partícula = st.selectbox("Selecione a partícula: ", ['Elétron', 'Próton', 'Muon', ])

potencial = st.selectbox("Selecione o potencial: ", ['Partícula Livre'])

if tempo == 'Dependente':
    if dimensão == '1D':
        x, t = symbols('x t')

        col1, col2, col3 = st.columns(3)

        with col1:
            Largura = st.number_input("Digite a Largura: ", value=8)

        with col2:
            n = st.number_input("Digite o número quântico: ", value=4)

        with col3:
            t_p = st.text_input("Digite o tempo (ex: 1e-20):", "1e-20")

        hbar = 1.054718e-34

        k = n * pi / Largura

        if partícula == 'Elétron':
            m = 9.1093837e-31

        elif partícula == 'Próton':
            m = 1.672621777e-27

        elif partícula == 'Muon':
            m = 1.9e-28

        if dimensão == '1D':
            if potencial == 'Partícula Livre':
                V = 0
                E = n**2*hbar**2*pi**2 / (2*m*Largura**2)
                X = exp(1j*k*x)

            elif potencial == 'Oscilador Hârmonico':
                from sympy import hermite
                o = st.number_input("Digite a frequência angular (sla oq é, mas digita aí que dá certo): ", value=1)
                V = 1/2 * m * o ** 2 * x ** 2
                E = hbar*o*(n+1/2)
                alpha = (sqrt(m*o/hbar)*x)
                X = hermite(n,alpha)*exp(-alpha**2/2)

        if st.button("Calcule"):
            try:
                C = sqrt(2/Largura)

                T = exp(-1j*E/hbar*t)
                
                PSI = C * X * T

                H_psi = ((-hbar**2 / (2*m)) * diff(PSI, x, 2) + V * PSI)

                ih_t = 1j*hbar*diff(PSI, t)

                st.title('Função de onda')
                st.write(PSI)

                psi_function = lambdify((x,t),PSI,"numpy")

                x_vals = np.linspace(0,float(Largura),10000)

                t_vals = float(t_p)

                z =  psi_function(x_vals, t_vals)

                fig, ax = plt.subplots()
                plt.style.use('dark_background')
                ax.plot(x_vals, np.real(z), label='Função de Onda (dependente do tempo), real', color='blue')

                ax.plot(x_vals, np.imag(z), label='Função de Onda (dependente do tempo), imag', color='red')

                ax.grid()

                ax.axhline(0,color="black")

                ax.axvline(0,color="black")

                ax.legend()

                ax.set_ylabel('RE(PSI) e IM(PSI)')

                ax.set_xlabel('Posição')

                st.pyplot(fig)
                    
            except Exception as e:
                st.error(f"ERROR: {e}")

    if dimensão == '2D':
        x, y, t = symbols('x y t')

        col1, col2, col3, col4, col5 = st.columns(5)

        with col1:
            Largura = st.number_input("Digite a Largura: ", value=1)

        with col2:
            Altura = st.number_input('Digite a Altura', value=1)


        with col3:
            nx = st.number_input("Digite o número quântico de x: ", value=1)

        with col4:
            ny = st.number_input("Digite o número quântico de y: ", value=1)

        with col5:
            t_p = st.number_input("Digite o instante de tempo: ", value=10.0, step=0.1)

        # H*psi(x,t) = 1j*h*d(psi)/dt
        # [-h²/2m * X''/X + V] = 1j*h*T'/T

        hbar = 1.054718e-34

        kx = nx * pi / Largura

        ky = ny * pi / Altura

        k = sqrt(kx**2+ky**2)

        if partícula == 'Elétron':
            m = 9.1093837e-31

        elif partícula == 'Próton':
            m = 1.672621777e-27

        elif partícula == 'Muon':
            m = 1.9e-28

        if dimensão == '2D':
            if potencial == 'Partícula Livre':
                V = 0
                E = k**2*hbar**2/(2*m)
                X = exp(1j*kx*x)
                Y = exp(1j*ky*y)

            elif potencial == 'Oscilador Hârmonico':
                from sympy import hermite
                o = st.number_input("Digite a frequência angular (sla oq é, mas digita aí que dá certo): ", value=1)
                V = 1/2 * m * o ** 2 * x ** 2 + 1/2 * m * o ** 2 * y ** 2
                E = hbar*o*(nx+ny+1)
                alpha = m*o/hbar
                X = hermite(nx, sqrt(alpha)*x) * exp( -alpha*x**2/2 )

                Y = hermite(ny, sqrt(alpha)*y) * exp( -alpha*y**2/2 )

        if st.button("Calcule"):
            try:
                C = sqrt(4/(Altura*Largura))

                T = exp( -1j * E * t)
                
                PSI = C * X * Y * T

                H_psi = ((-hbar**2 / (2*m)) * diff(PSI, x, 2)+diff(PSI, y, 2) + V * PSI)

                ih_t = 1j*hbar*diff(PSI, t)

                st.title('Função de onda')
                st.write(PSI)

                psi_function = lambdify((x,y,t),PSI,"numpy")

                x_vals = np.linspace(0,float(Largura),100)
                y_vals = np.linspace(0,float(Altura), 100)
                t_vals = t_p
                XX, YY = np.meshgrid(x_vals,y_vals)
                z =  psi_function(XX, YY, t_vals)
                
                fig = plt.figure(figsize=(10,7))

                ax = fig.add_subplot(111, projection='3d')

                ax.plot_surface(
                    XX, YY, np.real(z), 
                    label='Função de Onda, real', 
                    cmap='viridis',
                    edgecolor='none',
                    antialiased=True,
                    shade=True
                    )
                ax.set_xlabel('posição x')
                ax.set_ylabel('posição y')
                ax.set_zlabel('Função de Onda real')
                st.pyplot(fig)
                    
            except Exception as e:
                st.error(f"ERROR: {e}")

    if dimensão == '3D':
        x, y, z, t = symbols('x y z t')

        col1, col2, col3, col4, col5, col6 = st.columns(6)

        with col1:
            Largura = st.number_input("Digite a Largura: ", value=1)

        with col2:
            Altura = st.number_input('Digite a Altura', value=1)

        with col3:
            Profundidade = st.number_input("Digite a Profundidade", value=1)
        with col4:
            nx = st.number_input("Digite o número quântico de x: ", value=1)

        with col5:
            ny = st.number_input("Digite o número quântico de y: ", value=1)

        with col6:
            nz = st.number_input("Digite o número quantico de z: ", value=1)

        # H*psi(x,t) = 1j*h*d(psi)/dt
        # [-h²/2m * X''/X + V] = 1j*h*T'/T
        t_p = st.number_input("Digite o instante de tempo: ", value=10.0, step=0.1)

        hbar = 1.054718e-34

        kx = nx * pi / Largura

        ky = ny * pi / Altura

        kz = nz * pi / Profundidade

        k = sqrt(kx**2+ky**2+kz**2)

        if partícula == 'Elétron':
            m = 9.1093837e-31

        elif partícula == 'Próton':
            m = 1.672621777e-27

        elif partícula == 'Muon':
            m = 1.9e-28

        if dimensão == '3D':
            if potencial == 'Partícula Livre':
                V = 0
                E = k**2*hbar**2/(2*m)
                X = exp(1j*kx*x)
                Y = exp(1j*ky*y)
                Z = exp(1j*kz*z)

        if st.button("Calcule"):

            try:
                C = sqrt(8/(Largura*Altura*Profundidade))

                T = exp( -1j * E * t/hbar)
                
                PSI = C * X * Y * Z * T

                H_psi = ((-hbar**2 / (2*m)) * (diff(PSI, x, 2) + diff(PSI, y, 2)+ diff(PSI, z, 2)) + V * PSI)

                ih_t = 1j*hbar*diff(PSI, t)

                st.title('Equação de Schrondinger')
                H_psi=simplify(H_psi)
                ih_t=simplify(ih_t)
                st.write(H_psi, "=", ih_t) 

                st.title("Função de Onda")

                st.write(PSI)

            except Exception as e:
                st.error(f"Erro: {e}")

if tempo == 'Independente':
    if dimensão == '1D':
        x = symbols('x')

        col1, col2, col3 = st.columns(3)

        with col1:
            Largura = st.number_input("Digite a Largura: ", value=8)

        with col2:
            n = st.number_input("Digite o número quântico: ", value=4)

        hbar = 1.054718e-34

        k = n * pi / Largura

        if partícula == 'Elétron':
            m = 9.1093837e-31

        elif partícula == 'Próton':
            m = 1.672621777e-27

        elif partícula == 'Muon':
            m = 1.9e-28

        if dimensão == '1D':
            if potencial == 'Partícula Livre':
                V = 0
                E = n**2*hbar**2*pi**2 / (2*m*Largura**2)
                X = exp(1j*k*x)

            elif potencial == 'Oscilador Hârmonico':
                from sympy import hermite
                o = st.number_input("Digite a frequência angular (sla oq é, mas digita aí que dá certo): ", value=1)
                V = 1/2 * m * o ** 2 * x ** 2
                E = hbar*o*(n+1/2)
                alpha = (sqrt(m*o/hbar)*x)
                X = hermite(n,alpha)*exp(-alpha**2/2)

        if st.button("Calcule"):
            try:
                C = sqrt(2/Largura)

                PSI = C * X 

                H_psi = ((-hbar**2 / (2*m)) * diff(PSI, x, 2) + V * PSI)

                E_psi = E*PSI

                st.header(f"Equação de Schrodinger Independente do Tempo de {potencial}")
                st.write(H_psi, '=', E_psi)

                st.header('Função de onda independente do tempo')
                st.write(PSI)

                psi_function = lambdify(x,PSI,"numpy")

                x_vals = np.linspace(0,float(Largura),10000)
                y = psi_function(x_vals)


                fig, ax = plt.subplots()

                ax.plot(x_vals, np.real(y), label='Função de Onda, real', color='blue')

                ax.plot(x_vals, np.imag(y), label='Função de Onda, imag', color='red')

                ax.grid()

                ax.axhline(0,color="black")

                ax.axvline(0,color="black")

                ax.legend()

                ax.set_ylabel('RE(PSI) e IM(PSI)')

                ax.set_xlabel('Posição')

                st.pyplot(fig)
                    
            except Exception as e:
                st.error(f"ERROR: {e}")

    if dimensão == '2D':
        x, y = symbols('x y')

        col1, col2, col3, col4, col5 = st.columns(5)

        with col1:
            Largura = st.number_input("Digite a Largura: ", value=1)

        with col2:
            Altura = st.number_input('Digite a Altura', value=1)


        with col3:
            nx = st.number_input("Digite o número quântico de x: ", value=1)

        with col4:
            ny = st.number_input("Digite o número quântico de y: ", value=1)

        hbar = 1.054718e-34

        kx = nx * pi / Largura

        ky = ny * pi / Altura

        k = sqrt(kx**2+ky**2)

        if partícula == 'Elétron':
            m = 9.1093837e-31

        elif partícula == 'Próton':
            m = 1.672621777e-27

        elif partícula == 'Muon':
            m = 1.9e-28

        if dimensão == '2D':
            if potencial == 'Partícula Livre':
                V = 0
                E = k**2*hbar**2/(2*m)
                X = exp(1j*kx*x)
                Y = exp(1j*ky*y)

            elif potencial == 'Oscilador Hârmonico':
                from sympy import hermite
                o = st.number_input("Digite a frequência angular (sla oq é, mas digita aí que dá certo): ", value=1)
                V = 1/2 * m * o ** 2 * x ** 2 + 1/2 * m * o ** 2 * y ** 2
                E = hbar*o*(nx+ny+1)
                alpha = m*o/hbar
                X = hermite(nx, sqrt(alpha)*x) * exp( -alpha*x**2/2 )

                Y = hermite(ny, sqrt(alpha)*y) * exp( -alpha*y**2/2 )

        if st.button("Calcule"):
            try:
                C = sqrt(4/(Altura*Largura))

                PSI = C * X * Y

                H_psi = ((-hbar**2 / (2*m)) * diff(PSI, x, 2)+diff(PSI, y, 2) + V * PSI)

                E_psi = E * PSI

                st.header('Eq de Schrondigner Independente do tempo')
                st.write(H_psi, "=", E_psi)

                st.title('Função de onda')
                st.write(PSI)

                psi_function = lambdify((x,y),PSI,"numpy")

                x_vals = np.linspace(0,float(Largura),100)
                y_vals = np.linspace(0,float(Altura), 100)
                XX, YY = np.meshgrid(x_vals,y_vals)
                z =  psi_function(XX, YY)

                fig = plt.figure(figsize=(10,7))
                ax = fig.add_subplot(111, projection='3d')

                ax.plot_surface(
                    XX, YY, np.real(z), 
                    label='Função de Onda, real', 
                    cmap='viridis',
                    edgecolor='none',
                    antialiased=True,
                    shade=True
                    )
                ax.set_xlabel('posição x')
                ax.set_ylabel('posição y')
                ax.set_zlabel('Função de Onda real')
                st.pyplot(fig)
                    
            except Exception as e:
                st.error(f"ERROR: {e}")
    if dimensão == '3D':
        x, y, z = symbols('x y z')

        col1, col2, col3, col4, col5, col6 = st.columns(6)

        with col1:
            Largura = st.number_input("Digite a Largura: ", value=1)

        with col2:
            Altura = st.number_input('Digite a Altura', value=1)

        with col3:
            Profundidade = st.number_input("Digite a Profundidade", value=1)
        with col4:
            nx = st.number_input("Digite o número quântico de x: ", value=1)

        with col5:
            ny = st.number_input("Digite o número quântico de y: ", value=1)

        with col6:
            nz = st.number_input("Digite o número quantico de z: ", value=1)

        # H*psi(x,t) = 1j*h*d(psi)/dt
        # [-h²/2m * X''/X + V] = 1j*h*T'/T

        hbar = 1.054718e-34

        kx = nx * pi / Largura

        ky = ny * pi / Altura

        kz = nz * pi / Profundidade

        k = sqrt(kx**2+ky**2+kz**2)

        if partícula == 'Elétron':
            m = 9.1093837e-31

        elif partícula == 'Próton':
            m = 1.672621777e-27

        elif partícula == 'Muon':
            m = 1.9e-28

        if dimensão == '3D':
            if potencial == 'Partícula Livre':
                V = 0
                E = k**2*hbar**2/(2*m)
                X = exp(1j*kx*x)
                Y = exp(1j*ky*y)
                Z = exp(1j*kz*z)

        if st.button("Calcule"):

            try:
                C = sqrt(8/(Largura*Altura*Profundidade))

                PSI = C * X * Y * Z

                H_psi = ((-hbar**2 / (2*m)) * (diff(PSI, x, 2) + diff(PSI, y, 2)+ diff(PSI, z, 2)) + V * PSI)

                E_psi = E*PSI

                st.title('Equação de Schrondinger')
                H_psi=simplify(H_psi)
                E_psi=simplify(E_psi)
                st.write(H_psi, "=", E_psi) 

                st.title("Função de Onda")

                st.write(PSI)

            except Exception as e:
                st.error(f"Erro: {e}")

