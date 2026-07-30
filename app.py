import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import math

# ==========================================
# CONFIGURAÇÃO DA PÁGINA STREAMLIT
# ==========================================
st.set_page_config(page_title="Calculadora de Empuxo - Muro de Arrimo", layout="wide")

st.title("🧱 Calculadora de Empuxo em Muros de Arrimo")
st.markdown('''
Esta ferramenta calcula as pressões atuantes em um muro de arrimo utilizando a **Teoria de Rankine** 
e o **Princípio das Tensões Efetivas de Terzaghi** (considerando a presença de lençol freático).
''')

# ==========================================
# BARRA LATERAL (INPUTS DO USUÁRIO)
# ==========================================
st.sidebar.header("Parâmetros de Entrada")

# --- Geometria ---
st.sidebar.subheader("Geometria")
H = st.sidebar.number_input("Altura do Muro (H) em metros", min_value=1.0, max_value=20.0, value=4.0, step=0.5)

sem_agua = st.sidebar.checkbox("Lençol freático ausente (Muro Drenado)", value=True)
if sem_agua:
    zw = float(H)
else:
    zw = st.sidebar.number_input("Profundidade do Lençol Freático (m)", min_value=0.0, max_value=float(H), value=2.0, step=0.5, help="Medido a partir do topo do muro.")

# --- Estrutura do Muro ---
st.sidebar.subheader("Estrutura do Muro")
tipo_muro = st.sidebar.radio(
    "Tipo de Lançamento de Cargas", 
    ["Maciço / Contínuo (Carga linear em kN/m)", "Reticulado com Pilares (Carga pontual em kN)"]
)

if tipo_muro == "Reticulado com Pilares (Carga pontual em kN)":
    dist_pilares = st.sidebar.number_input("Distância entre Pilares (l) em metros", min_value=0.5, max_value=10.0, value=2.0, step=0.1)
else:
    dist_pilares = 1.0  # Para muro contínuo, analisamos uma faixa de 1 metro linear

# --- Propriedades do Solo ---
st.sidebar.subheader("Propriedades do Solo")
tipo_solo = st.sidebar.radio(
    "Comportamento do Solo",
    ["Não Coesivo (Areias, Cascalhos)", "Coesivo (Argilas, Siltes)"]
)

gamma = st.sidebar.number_input("Peso Específico Natural (kN/m³)", min_value=10.0, value=16.0, step=0.5)

if sem_agua:
    gamma_sat = gamma
else:
    st.sidebar.caption("Água presente: informe o peso do solo saturado.")
    gamma_sat = st.sidebar.number_input("Peso Específico Saturado (kN/m³)", min_value=10.0, value=19.0, step=0.5)

phi = st.sidebar.number_input("Ângulo de Atrito Interno (graus)", min_value=1.0, value=30.0, step=1.0)

if tipo_solo == "Coesivo (Argilas, Siltes)":
    c = st.sidebar.number_input("Coesão (kPa)", min_value=0.0, value=10.0, step=1.0)
else:
    c = 0.0

# --- Condições Adicionais ---
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
if zw == 0 or zw == H:
    z_pts = np.array([0.0, H])
else:
    z_pts = np.array([0.0, zw, H])

sigma_v_efetiva = np.zeros_like(z_pts)
u_agua = np.zeros_like(z_pts)

# Cálculo das tensões ponto a ponto
for i, z in enumerate(z_pts):
    if z <= zw:
        sigma_v_efetiva[i] = q + (gamma * z)
        u_agua[i] = 0.0
    else:
        gamma_submerso = gamma_sat - gamma_w
        sigma_v_efetiva[i] = q + (gamma * zw) + (gamma_submerso * (z - zw))
        u_agua[i] = gamma_w * (z - zw)

# Pressões horizontais
pressao_solo_efetiva = (sigma_v_efetiva * ka) - (2 * c * math.sqrt(ka))
pressao_solo_efetiva = np.maximum(pressao_solo_efetiva, 0) # Ignora zona de tração
pressao_total = pressao_solo_efetiva + u_agua

# Integração (Cálculo da força em kN/m)
empuxo_solo = np.trapezoid(pressao_solo_efetiva, z_pts)
empuxo_agua = np.trapezoid(u_agua, z_pts)
empuxo_total = empuxo_solo + empuxo_agua

# Cálculo do Ponto de Aplicação (A partir da base)
if empuxo_total > 0:
    alturas_da_base = H - z_pts
    momento_estatico = np.trapezoid(pressao_total * alturas_da_base, z_pts)
    y_aplicacao_base = momento_estatico / empuxo_total
else:
    y_aplicacao_base = 0.0

# ==========================================
# EXIBIÇÃO DOS RESULTADOS (UI STREAMLIT)
# ==========================================
st.divider()
st.subheader("Resultados do Empuxo")

col1, col2, col3 = st.columns(3)
col1.metric("Ponto de Aplicação (da base)", f"{y_aplicacao_base:.2f} m")
col2.metric("Coeficiente Ativo (Ka)", f"{ka:.3f}")
col3.metric("Empuxo Total Distribuído", f"{empuxo_total:.1f} kN/m")

# Se for muro reticulado, damos ênfase na carga pontual final por pilar
if tipo_muro == "Reticulado com Pilares (Carga pontual em kN)":
    st.info(f"💡 **Carga Pontual no Pilar:** Como a distância entre pilares é de {dist_pilares}m, a força horizontal resultante a ser aplicada em cada pilar no software estrutural (ex: Eberick) é de **{(empuxo_total * dist_pilares):.1f} kN**, localizada a **{y_aplicacao_base:.2f} m** da base.")

# Gráfico
st.subheader("Diagrama de Pressões Horizontais")
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(pressao_solo_efetiva, z_pts, color='saddlebrown', linewidth=2, label='Pressão do Solo (Efetiva)')
ax.plot(pressao_total, z_pts, color='black', linewidth=2, linestyle='--', label='Pressão Total')
ax.fill_betweenx(z_pts, 0, pressao_solo_efetiva, color='saddlebrown', alpha=0.3)
ax.fill_betweenx(z_pts, pressao_solo_efetiva, pressao_total, color='dodgerblue', alpha=0.5, label='Pressão da Água (Hidrostática)')
if 0 < zw < H:
    ax.axhline(y=zw, color='blue', linestyle='-.', linewidth=1.5, alpha=0.7)
    ax.text(pressao_total.max() * 0.8, zw - 0.1, f'NA (z={zw}m)', color='blue', fontsize=10)
ax.invert_yaxis()
ax.set_ylabel('Profundidade (m)')
ax.set_xlabel('Pressão Horizontal (kPa)')
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='lower left')
st.pyplot(fig)
