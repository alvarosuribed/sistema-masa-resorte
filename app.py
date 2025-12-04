import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Configuración de la página
st.set_page_config(page_title="Sistema Masa-Resorte-Amortiguador", layout="wide")

# CSS para eliminar espacio superior y ajustar padding
st.markdown("""
    <style>
        .block-container {
            padding-top: 1rem;
            padding-bottom: 0rem;
        }
        h1 {
            margin-top: 0rem;
        }
        .stLatex {
            font-size: 1.1rem;
        }
    </style>
    """, unsafe_allow_html=True)

# Título
st.title("Sistema Masa-Resorte-Amortiguador Interactivo")
st.markdown("**Problema 18:** Análisis dinámico de varilla articulada con forzamiento sinusoidal")
st.caption("Simulador que resuelve numéricamente la ecuación diferencial y permite validar el modelo teórico mediante análisis paramétrico en tiempo real.")
st.markdown("---")

# Sidebar con parámetros
st.sidebar.header("⚙️ Parámetros del Sistema")

m = st.sidebar.slider("Masa m (kg)", 0.5, 10.0, 2.0, 0.5)
k = st.sidebar.slider("Rigidez k (N/m)", 10.0, 500.0, 100.0, 10.0)
c = st.sidebar.slider("Amortiguamiento c (N·s/m)", 0.1, 100.0, 10.0, 1.0)
F = st.sidebar.slider("Amplitud Fuerza F (N)", 1.0, 50.0, 10.0, 1.0)
omega = st.sidebar.slider("Frecuencia ω (rad/s)", 0.1, 20.0, 5.0, 0.1)

t_max = st.sidebar.slider("Tiempo de simulación (s)", 5.0, 30.0, 15.0, 1.0)

# ==========================================
# EDO DINÁMICA CON VALORES ACTUALES
# ==========================================
st.markdown("### 📐 Ecuación Diferencial del Sistema")

# Función para formatear números de forma limpia
def formato_num(valor, decimales=1):
    """Formatea número eliminando decimales innecesarios."""
    if valor == int(valor):
        return str(int(valor))
    return f"{valor:.{decimales}f}"

# Calcular coeficientes
coef_aceleracion = 4 * m
coef_velocidad = c
coef_posicion = k
coef_fuerza = 4 * F

# Construir la ecuación LaTeX dinámica
edo_latex = rf"{formato_num(coef_aceleracion)}\ddot{{x}} + {formato_num(coef_velocidad)}\dot{{x}} + {formato_num(coef_posicion)}x = {formato_num(coef_fuerza)}\sin({formato_num(omega)}t)"

# Mostrar ecuación general y con valores
col_edo1, col_edo2 = st.columns(2)

with col_edo1:
    st.markdown("**Forma general:**")
    st.latex(r"4m\ddot{x} + c\dot{x} + kx = 4F\sin(\omega t)")

with col_edo2:
    st.markdown("**Con valores actuales:**")
    st.latex(edo_latex)

# Mostrar forma normalizada (dividida entre el coeficiente de ẍ)
st.markdown(r"**Forma normalizada** (dividida entre el coeficiente de $\ddot{x}$):")
a1 = coef_velocidad / coef_aceleracion
a0 = coef_posicion / coef_aceleracion
f0 = coef_fuerza / coef_aceleracion

edo_normalizada = rf"\ddot{{x}} + {formato_num(a1, 2)}\dot{{x}} + {formato_num(a0, 2)}x = {formato_num(f0, 2)}\sin({formato_num(omega)}t)"
st.latex(edo_normalizada)

st.markdown("---")

# ==========================================
# CÁLCULO DE PARÁMETROS DEL SISTEMA
# ==========================================
omega_0 = 0.5 * np.sqrt(k / m)
zeta = c / (4 * np.sqrt(m * k))
c_critico = 4 * np.sqrt(k * m)

# Cálculo de omega_d
if zeta < 1:  # Subamortiguado
    omega_d = omega_0 * np.sqrt(1 - zeta**2)
else:
    omega_d = 0

# Clasificación de amortiguamiento
if zeta < 1:
    tipo_amort, color_amort = "Subamortiguado", "🟢"
elif abs(zeta - 1) < 0.01:
    tipo_amort, color_amort = "Críticamente amortiguado", "🟡"
else:
    tipo_amort, color_amort = "Sobreamortiguado", "🔴"

# ==========================================
# MOSTRAR FRECUENCIAS DINÁMICAMENTE
# ==========================================
st.markdown("### 📊 Frecuencias del Sistema")

col_freq1, col_freq2, col_freq3 = st.columns(3)

with col_freq1:
    st.latex(rf"\omega_0 = \frac{{1}}{{2}}\sqrt{{\frac{{k}}{{m}}}} = \frac{{1}}{{2}}\sqrt{{\frac{{{formato_num(k)}}}{{{formato_num(m)}}}}} = {omega_0:.3f} \text{{ rad/s}}")

with col_freq2:
    if omega_d > 0:
        st.latex(rf"\omega_d = \omega_0\sqrt{{1-\zeta^2}} = {omega_0:.3f}\sqrt{{1-{zeta:.3f}^2}} = {omega_d:.3f} \text{{ rad/s}}")
    else:
        st.latex(r"\omega_d = 0 \text{ (No oscila)}")

with col_freq3:
    st.latex(rf"\zeta = \frac{{c}}{{4\sqrt{{mk}}}} = \frac{{{formato_num(c)}}}{{4\sqrt{{{formato_num(m)} \cdot {formato_num(k)}}}}} = {zeta:.4f}")

# Métricas compactas
st.markdown("---")
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("ω₀ (Natural)", f"{omega_0:.2f} rad/s")
with col2:
    st.metric("ωd (Amortiguada)", f"{omega_d:.2f} rad/s" if omega_d > 0 else "No oscila")
with col3:
    st.metric("ζ (Razón amort.)", f"{zeta:.3f}")
with col4:
    st.metric("Tipo", f"{color_amort} {tipo_amort}")

# Advertencias
if omega_d > 0 and abs(omega - omega_d) < 0.5:
    st.warning("⚠️ **ADVERTENCIA:** Frecuencia de excitación cerca de resonancia! Amplitudes grandes.")

st.info(f"ℹ️ **Amortiguamiento crítico:** c_crítico = 4√(mk) = 4√({formato_num(m)}·{formato_num(k)}) = {c_critico:.2f} N·s/m")

# ==========================================
# RESOLVER EDO
# ==========================================
@st.cache_data
def resolver_edo(m, c, k, F, omega, t_max):
    def sistema_edo(y, t):
        x, v = y
        dxdt = v
        dvdt = (4 * F * np.sin(omega * t) - c * v - k * x) / (4 * m)
        return [dxdt, dvdt]
    
    t = np.linspace(0, t_max, 2000)
    sol = odeint(sistema_edo, [0.0, 0.0], t)
    return t, sol[:, 0], sol[:, 1]

t, x, v = resolver_edo(m, c, k, F, omega, t_max)

# ==========================================
# GRÁFICAS
# ==========================================
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 7))

# Gráfica 1: Desplazamiento
ax1.plot(t, x, 'b-', linewidth=2, label='x(t)')
ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax1.set_xlabel('Tiempo (s)', fontsize=10)
ax1.set_ylabel('Desplazamiento x(t) (m)', fontsize=10)
ax1.set_title('Desplazamiento del extremo C', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=9)

# Gráfica 2: Velocidad
ax2.plot(t, v, 'r-', linewidth=2, label='v(t) = dx/dt')
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax2.set_xlabel('Tiempo (s)', fontsize=10)
ax2.set_ylabel('Velocidad v(t) (m/s)', fontsize=10)
ax2.set_title('Velocidad del extremo C', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=9)

# Gráfica 3: Diagrama de fase
ax3.plot(x, v, 'g-', linewidth=1.5, alpha=0.7)
ax3.plot(x[0], v[0], 'go', markersize=8, label='Inicio')
ax3.plot(x[-1], v[-1], 'ro', markersize=8, label='Final')
ax3.set_xlabel('Desplazamiento x (m)', fontsize=10)
ax3.set_ylabel('Velocidad v (m/s)', fontsize=10)
ax3.set_title('Diagrama de Fase', fontsize=12, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.legend(fontsize=9)
ax3.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax3.axvline(x=0, color='k', linestyle='--', alpha=0.3)

plt.tight_layout()
st.pyplot(fig)

# ==========================================
# ANÁLISIS ADICIONAL
# ==========================================
st.markdown("---")
st.markdown("### 📈 Análisis de Resultados")

col1, col2 = st.columns(2)

with col1:
    st.markdown("**Parámetros actuales:**")
    st.markdown(f"""
    | Parámetro | Valor |
    |-----------|-------|
    | Masa efectiva (4m) | {coef_aceleracion:.2f} kg |
    | Rigidez (k) | {k:.2f} N/m |
    | Amortiguamiento (c) | {c:.2f} N·s/m |
    | c/c_crítico | {c/c_critico:.3f} |
    | Fuerza (4F) | {coef_fuerza:.2f} N |
    | Frecuencia (ω) | {omega:.2f} rad/s |
    """)

with col2:
    amplitud_max = np.max(np.abs(x[int(len(x)*0.5):]))
    periodo = 2*np.pi/omega_d if omega_d > 0 else float('inf')
    relacion_omega = f"{omega/omega_d:.3f}" if omega_d > 0 else "N/A"
    
    st.markdown("**Resultados de la simulación:**")
    st.markdown(f"""
    | Característica | Valor |
    |----------------|-------|
    | Amplitud máx. (estado estable) | {amplitud_max:.4f} m |
    | Período de oscilación | {periodo:.2f} s |
    | Relación ω/ωd | {relacion_omega} |
    | Frecuencia excitación | {omega:.2f} rad/s |
    | Frecuencia natural | {omega_0:.2f} rad/s |
    """)

# Footer
st.markdown("---")
st.markdown("**💡 Instrucciones:** Ajusta los sliders en la barra lateral para ver cómo cambian la ecuación y las frecuencias en tiempo real.")
st.markdown(f"**🎯 Tip:** Para ver sobreamortiguamiento, aumenta c por encima de {c_critico:.1f} N·s/m")