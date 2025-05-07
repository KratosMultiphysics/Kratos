import numpy as np
import plotly.graph_objects as go


def plot_sigma(sigma_1, sigma_3):
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=sigma_3,
        y=sigma_1,
        mode='markers',
        marker=dict(size=10, color='blue'),
        name='σ₁ vs σ₃'
    ))

    fig.update_layout(
        title=dict(
            text='σ₁ vs σ₃ Plot',
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title='σ₃ (Principal Stress 3) [kN/m²]',
            showline=True,
            autorange='reversed',
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        yaxis=dict(
            title=' σ₁ (Principal Stress 1) [kN/m²]',
            showline=True,
            autorange='reversed',
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        template='plotly_white',
    )
    fig.update_layout(
        xaxis=dict(rangemode='tozero'),
        yaxis=dict(rangemode='tozero'),
    )
    return fig


def plot_delta_sigma(vertical_strain, sigma_diff):
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=vertical_strain,
        y=sigma_diff,
        mode='markers',
        marker=dict(size=10, color='blue'),
        name='|σ₁ - σ₃|'
    ))

    fig.update_layout(
        title=dict(
            text='|σ₁ - σ₃| vs Vertical Strain Plot',
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title='Vertical Strain [-]',
            autorange='reversed',
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        yaxis=dict(
            title='|σ₁ - σ₃| [kN/m²]',
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        template='plotly_white'
    )

    fig.update_layout(
        xaxis=dict(rangemode='tozero'),
        yaxis=dict(rangemode='tozero')
    )
    return fig

def plot_volumetric_strain(vertical_strain, volumetric_strain):
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=vertical_strain,
        y=volumetric_strain,
        mode='markers',
        marker=dict(size=10, color='blue'),
        name='Volumetric Strain'
    ))

    fig.update_layout(
        title=dict(
            text='Volumetric Strain vs Vertical Strain Plot',
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title='Vertical Strain [−]',
            autorange='reversed',
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        yaxis=dict(
            title='Volumetric Strain [−]',
            autorange='reversed',
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        template='plotly_white'
    )

    fig.update_layout(
        xaxis=dict(rangemode='tozero'),
        yaxis=dict(rangemode='tozero')
    )
    return fig

def plot_mohr_coulomb_circle(sigma_1, sigma_3, cohesion, friction_angle):
    center = (sigma_1 + sigma_3) / 2
    radius = (sigma_1 - sigma_3) / 2
    theta = np.linspace(0, np.pi, 200)
    sigma = center + radius * np.cos(theta)
    tau = radius * np.sin(theta)

    phi_rad = np.radians(friction_angle)
    x_line = np.linspace(0, sigma_1, 200)
    y_line = x_line * np.tan(phi_rad) - cohesion

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=sigma,
        y=tau,
        mode='lines',
        name='Mohr-Coulomb Circle',
        line=dict(color='blue', width=2)
    ))

    fig.add_trace(go.Scatter(
        x=x_line,
        y=y_line,
        mode='lines',
        name="Mobilized Shear Stress = σ' * tan(ϕ°) + c",
        line=dict(color='red', width=2, dash='dash')
    ))

    fig.update_layout(
        title=dict(
            text="Mohr-Coulomb Circle",
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title="σ' (Effective Stress) [kN/m²]",
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True,
            autorange='reversed'
        ),
        yaxis=dict(
            title="τ (Mobilized Shear Stress) [kN/m²]",
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True,
            autorange='reversed'
        ),
        template='plotly_white',
        xaxis_scaleanchor="y"
    )
    return fig

def plot_p_q(p_list, q_list):
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=p_list,
        y=q_list,
        mode='markers+lines',
        name="p' vs q",
        marker=dict(size=8, color='blue')
    ))

    fig.update_layout(
        title=dict(text="Mean Effective Stress vs Deviatoric Stress", x=0.5),
        xaxis=dict(
            title="p' (Mean Effective Stress) [kN/m²]",
            autorange='reversed',
            showline=True,
            mirror=True,
            linecolor='black'
        ),
        yaxis=dict(
            title="q (Deviatoric Stress) [kN/m²]",
            showline=True,
            mirror=True,
            linecolor='black'
        ),
        template='plotly_white'
    )

    fig.update_layout(
        xaxis=dict(rangemode='tozero'),
        yaxis=dict(rangemode='tozero')
    )
    return fig