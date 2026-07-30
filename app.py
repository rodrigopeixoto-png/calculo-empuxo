import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import math

# ==========================================
# CONFIGURAÇÃO DA PÁGINA STREAMLIT
# ==========================================
st.set_page_config(page_title="Calculadora de Empuxo - Muro de Arrimo", layout="wide")

st.title("🧱 Calculadora de Empuxo em Muros de Arrimo")
st.markdown("""
Esta ferramenta calcula as pressões atuantes em um muro de arrimo utilizando a **Teoria de Rankine** 
e o **Princípio das Tensões Efetivas de Terzaghi** (considerando a presença de lençol freático).
""")

# ==========================================
# BARRA LATERAL (INPUTS DO USUÁRIO)
# ==========================================
st.sidebar.header("Parâmetros de Entrada")

st.sidebar.subheader("Geometria")
H = st.sidebar.number_input("Altura do Muro (H) em metros", min_value=1.0, max_value=20.0, value=4.0, step=0.5)

# Checkbox para ausência de água
sem_agua = st.sidebar.checkbox("Lençol freático ausente (Muro Drenado)", value=False)

if sem_agua:
    zw = float(H)  # Se não há água, consideramos a profundidade igual à base (sem empuxo)
else:
    zw = st.sidebar.number_input("Profundidade do Lençol Freático (m)", min_value=0.0, max_value=float(H), value=2.0, step=0.5, help="Medido a partir do topo do muro.")

st.sidebar.subheader("Propriedades do Solo")
gamma = st.sidebar.number_input("Peso Específico Natural (kN/m³)", min_value=10.0, value=16.0, step=0.5)
gamma_sat = st.sidebar.number_input("Peso Específico Saturado (kN/m³)", min_value=10.0, value=19.0, step=0.5)
phi = st.sidebar.number_input("Ângulo de Atrito Interno (graus)", min_value=1.0, value=30.0, step=1.0)
c = st.sidebar.number_input("Coesão (kPa)", min_value=0.0, value=0.0, step=1.0)

st.sidebar.subheader("Condições Adicionais")
q = st.sidebar.number_input("Sobrecarga no Topo (kPa)", min_value=0.0, value=0.0, step=5.0)
beta = st.sidebar.number_input("Inclinação do Maciço - Beta (graus)", min_value=0.0, max_value=float(phi)-0.1, value=0.0, step=1.0)

# ==========================================
# LÓGICA DE CÁLCULO
# ==========================================
gamma_w = 9.81
beta_rad = math.radians(beta)
phi_rad = math.radians(phi)

# Coeficiente de Empuxo Ativo (Ka)
cos_beta = math.cos(beta_rad)
cos_phi = math.cos(phi_rad)
raiz = math.sqrt(max(0, cos_beta**2 - cos_phi**2))
ka = cos_beta * ((cos_beta - raiz) / (cos_beta + raiz))

# Arrays de profundidade (z)
# z=0 (topo), z=zw (nível da água), z=H (base)
if zw == 0:
    z_pts = np.array([0.0, H])
elif zw == H:
    z_pts = np.array([0.0, H])
else:
    z_pts = np.array([0.0, zw, H])

sigma_v_efetiva = np.zeros_like(z_pts)
u_agua = np.zeros_like(z_pts)

# Cálculo ponto a ponto
for i, z in enumerate(z_pts):
    if z <= zw:
        # Acima do nível d'água
        sigma_v_efetiva[i] = q + (gamma * z)
        u_agua[i] = 0.0
    else:
        # Abaixo do nível d'água
        gamma_submerso = gamma_sat - gamma_w
        sigma_v_efetiva[i] = q + (gamma * zw) + (gamma_submerso * (z - zw))
        u_agua[i] = gamma_w * (z - zw)

# Pressões horizontais
pressao_solo_efetiva = (sigma_v_efetiva * ka) - (2 * c * math.sqrt(ka))
pressao_solo_efetiva = np.maximum(pressao_solo_efetiva, 0) # Zona de tração não contribui com empuxo negativo

pressao_total = pressao_solo_efetiva + u_agua

# Integração (Cálculo da força total em kN/m) - Área do trapézio/triângulo
empuxo_solo = np.trapezoid(pressao_solo_efetiva, z_pts)
empuxo_agua = np.trapezoid(u_agua, z_pts)
empuxo_total = empuxo_solo + empuxo_agua

# ==========================================
# EXIBIÇÃO DOS RESULTADOS (UI STREAMLIT)
# ==========================================
col1, col2, col3, col4 = st.columns(4)
col1.metric("Coeficiente Ativo (Ka)", f"{ka:.3f}")
col2.metric("Empuxo Efetivo Solo", f"{empuxo_solo:.1f} kN/m")
col3.metric("Empuxo da Água", f"{empuxo_agua:.1f} kN/m")
col4.metric("Empuxo Total", f"{empuxo_total:.1f} kN/m")

st.divider()

# Gráfico
st.subheader("Diagrama de Pressões Horizontais")

fig, ax = plt.subplots(figsize=(10, 5))

# Plotagem
ax.plot(pressao_solo_efetiva, z_pts, color='saddlebrown', linewidth=2, label='Pressão do Solo (Efetiva)')
ax.plot(pressao_total, z_pts, color='black', linewidth=2, linestyle='--', label='Pressão Total')

# Preenchimento
ax.fill_betweenx(z_pts, 0, pressao_solo_efetiva, color='saddlebrown', alpha=0.3)
ax.fill_betweenx(z_pts, pressao_solo_efetiva, pressao_total, color='dodgerblue', alpha=0.5, label='Pressão da Água (Hidrostática)')

# Linha do Nível d'água
if 0 < zw < H:
    ax.axhline(y=zw, color='blue', linestyle='-.', linewidth=1.5, alpha=0.7)
    ax.text(pressao_total.max() * 0.8, zw - 0.1, f'NA (z={zw}m)', color='blue', fontsize=10)

ax.invert_yaxis()
ax.set_ylabel('Profundidade (m)')
ax.set_xlabel('Pressão Horizontal (kPa)')
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='lower left')

# Renderiza no Streamlit
st.pyplot(fig)

st.caption("Nota: Valores de empuxo negativo (zona de tração em solos coesivos) foram limitados a zero na plotagem e na integração de cálculo, conforme prática comum de projeto.")
