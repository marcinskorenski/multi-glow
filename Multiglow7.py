import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy.integrate import simpson
import numpy as np

st.set_page_config(layout="wide")

def calculate_auc(y_values, x_values):
    return simpson(y_values, x=x_values)

# Initialize session state
if 'analyses' not in st.session_state:
    st.session_state.analyses = []
if 'color_scale' not in st.session_state:
    st.session_state.color_scale = "Blues"
if 'show_main_app' not in st.session_state:
    st.session_state.show_main_app = False  # Start with welcome screen

# Welcome Window
def show_welcome_window():
    # Create two columns for title and image, adjusted ratio for larger image
    col1, col2 = st.columns([1.5, 1])  # Adjusted from [2, 1] to give image more space
    
    with col1:
        st.title("Welcome to Chemiluminescent Substrates Data Analysis")
    
    with col2:
        try:
            st.image("IMGP.jpg", width=300)  # Increased from 150 to 300
        except FileNotFoundError:
            st.error("Image file 'IMGP.jpg' not found. Please ensure it’s in the correct directory.")

    st.markdown("""
        ### About This Tool
        This platform allows users to explore and analyze data generated during the Multi Glow MSCA Fellowship project.  
        The data available here has already been presented at two conferences, 
        and we anticipate publishing our findings later this year. 
        
    """)
    if st.button("Let's Start", key="start_button"):
        st.session_state.show_main_app = True
        st.rerun()

# Main App Function
def show_main_app():
    # Sidebar configuration for enzyme selection with Multi-Select in a Form
    st.sidebar.title('Select Enzymes for Analysis')
    enzyme_options = ['TMPRSS2', 'uPA', 'HAT', 'FXa', 'Matriptase', 'Thrombin', 'Trypsin']

    with st.sidebar.form(key='enzyme_selection_form'):
        selected_enzymes = st.multiselect(
            'Choose Enzymes', 
            enzyme_options,
            default=st.session_state.analyses,
            key='enzyme_multiselect'
        )
        submit_button = st.form_submit_button(label="Update Analysis")

    if submit_button:
        if selected_enzymes != st.session_state.analyses:
            st.session_state.analyses = selected_enzymes
            st.rerun()

    st.sidebar.markdown("### Selected Analyses")
    to_remove = None
    for i, enzyme in enumerate(st.session_state.analyses):
        st.sidebar.write(f"{i+1}. {enzyme}")
        if st.sidebar.button(f"Remove {enzyme}", key=f'remove_{i}'):
            to_remove = i

    if st.sidebar.button("Reset All"):
        st.session_state.analyses = []
        st.rerun()

    if to_remove is not None:
        st.session_state.analyses.pop(to_remove)
        st.rerun()

    @st.cache_data
    def load_data(enzyme):
        file_path = f'{enzyme}.csv'
        with st.spinner(f"Loading data for {enzyme}..."):
            try:
                data = pd.read_csv(file_path, header=None, skiprows=4)
                time_points = data.iloc[:, 0].astype(float)
                substrates = {}
                columns_per_substrate = 3
                substrate_names = [f'Substrate {i+1}' for i in range(7)]
                substrate_sequences = {  
                    'Substrate 1': 'PEG-Ile-Gln-Phe-Arg-LUMI',
                    'Substrate 2': 'PEG-Gly-Thr-Ala-Arg-LUMI',
                    'Substrate 3': 'PEG-Arg-Gln-Asp-Arg-LUMI',
                    'Substrate 4': 'PEG-DPro-hArg-Nal-Arg-LUMI',
                    'Substrate 5': 'PEG-Arg-Gln-Arg-Arg-LUMI',
                    'Substrate 6': 'PEG-Nle-Lys-Pro-Arg-LUMI',
                    'Substrate 7': 'PEG-Gln-Ala-Arg-LUMI',
                }
                for i, substrate in enumerate(substrate_names):
                    start_col = 1 + i * columns_per_substrate
                    triplicates = data.iloc[:, start_col:start_col + columns_per_substrate]
                    substrates[substrate] = triplicates.median(axis=1)
                return time_points, substrates, substrate_sequences
            except FileNotFoundError:
                st.error(f"Data file for {enzyme} not found.")
                return None, None, None

    st.title("Chemiluminescent Data Analysis")

    # Initialize lists for heatmap data and substrate/enzyme tracking
    heatmap_data = []
    substrate_list = sorted([f'Substrate {i+1}' for i in range(7)])
    enzyme_list = st.session_state.analyses.copy()
    sequence_map = {}

    color_options = [
        "Blues", "Greens", "Reds", "Purples", "Oranges", "Viridis", "Plasma", "Inferno", "Magma", "Cividis",
        "RdBu", "PiYG", "PRGn", "BrBG", "PuOr", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Coolwarm", "Bwr"
    ]

    if not st.session_state.analyses:
        st.info("No data selected for analysis. Please select enzymes from the sidebar.")
    else:
        grid = st.columns(2) if len(st.session_state.analyses) > 1 else [st.container()]
        
        for idx, enzyme in enumerate(st.session_state.analyses):
            col = idx % 2
            with grid[col if len(st.session_state.analyses) > 1 else 0]:
                time_points, substrate_data, substrate_sequences = load_data(enzyme)
                
                if time_points is not None:
                    fig = go.Figure()
                    total_light_emitted = []
                    max_time = time_points.max()
                    animation_duration = 10000
                    num_frames = 50
                    frame_times = np.linspace(0, max_time, num_frames)

                    for substrate, y_values in substrate_data.items():
                        sequence = substrate_sequences.get(substrate, "Unknown Sequence")
                        hover_text = f"{substrate}: {sequence}"
                        fig.add_trace(go.Scatter(
                            x=time_points, 
                            y=y_values, 
                            mode='lines', 
                            name=substrate, 
                            hoverinfo='text', 
                            text=hover_text
                        ))
                        auc = calculate_auc(y_values, time_points)
                        total_light_emitted.append((substrate, auc, sequence))

                    frames = []
                    for t in frame_times:
                        frame_data = []
                        for substrate, y_values in substrate_data.items():
                            mask = time_points <= t
                            frame_x = time_points[mask]
                            frame_y = y_values[mask]
                            sequence = substrate_sequences.get(substrate, "Unknown Sequence")
                            hover_text = f"{substrate}: {sequence}"
                            frame_data.append(go.Scatter(
                                x=frame_x,
                                y=frame_y,
                                mode='lines',
                                name=substrate,
                                hoverinfo='text',
                                text=hover_text
                            ))
                        frames.append(go.Frame(data=frame_data))

                    fig.frames = frames
                    fig.update_layout(
                        title=f'Chemiluminescent Curves for {enzyme}',
                        xaxis_title='Time (seconds)',
                        yaxis_title='RLU',
                        width=1200 if len(st.session_state.analyses) == 1 else 600,
                        height=400,
                        updatemenus=[dict(
                            type="buttons",
                            direction="left",
                            buttons=[
                                dict(args=[None, {"frame": {"duration": animation_duration / num_frames, "redraw": True},
                                                "fromcurrent": True, "mode": "immediate"}],
                                     label="Play",
                                     method="animate"),
                                dict(args=[[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                                     label="Pause",
                                     method="animate")
                            ],
                            pad={"r": 10, "t": 10},
                            showactive=True,
                            x=0.1,
                            xanchor="right",
                            y=1.0,
                            yanchor="bottom",
                            font=dict(color="black"),
                            bgcolor="white"
                        )],
                        xaxis=dict(range=[0, max_time]),
                        yaxis=dict(range=[0, max(y_values.max() for y_values in substrate_data.values())])
                    )
                    st.plotly_chart(fig, use_container_width=True, key=f"curve_{enzyme}")

                    # Bar chart with percentage calculation
                    auc_df = pd.DataFrame(total_light_emitted, columns=['Substrate', 'Total Light Emitted', 'Sequence'])
                    max_auc = auc_df['Total Light Emitted'].max()
                    auc_df['Percentage'] = (auc_df['Total Light Emitted'] / max_auc) * 100
                    auc_df['Hover Text'] = auc_df.apply(lambda row: f"{row['Substrate']} ({row['Sequence']}): {row['Total Light Emitted']:.2f} ({row['Percentage']:.2f}%)", axis=1)

                    fig_bar = px.bar(
                        auc_df, 
                        x='Substrate', 
                        y='Total Light Emitted', 
                        color='Percentage',
                        color_continuous_scale=st.session_state.color_scale,
                        title=f'Total Light Emitted (AUC) for {enzyme}',
                        hover_data=['Hover Text']
                    )
                    fig_bar.update_traces(hovertemplate="%{customdata[0]}")
                    fig_bar.update_layout(coloraxis_showscale=True)
                    st.plotly_chart(fig_bar, use_container_width=True, key=f"bar_{enzyme}")
                    
                    # Collect data for heatmap
                    for _, row in auc_df.iterrows():
                        heatmap_data.append([row['Substrate'], enzyme, row['Percentage'], row['Sequence']])
                        sequence_map[row['Substrate']] = row['Sequence']

        # Color Theme Selection with Form
        if st.session_state.analyses:
            with st.form(key='color_form'):
                new_color_scale = st.selectbox(
                    "Select Color Theme", 
                    color_options, 
                    index=color_options.index(st.session_state.color_scale),
                    key="global_color_select"
                )
                submit_button_color = st.form_submit_button(label="Apply Color Theme")
            
            if submit_button_color:
                st.session_state.color_scale = new_color_scale
                st.rerun()

        # Generate and display a heatmap
        if heatmap_data:
            heatmap_df = pd.DataFrame(heatmap_data, columns=['Substrate', 'Enzyme', 'Percentage', 'Sequence'])
            
            # Ensure all substrates are present for each enzyme
            all_combinations = []
            for enzyme in enzyme_list:
                for substrate in substrate_list:
                    if not ((heatmap_df['Enzyme'] == enzyme) & (heatmap_df['Substrate'] == substrate)).any():
                        sequence = sequence_map.get(substrate, "Unknown Sequence")
                        all_combinations.append([substrate, enzyme, 0.0, sequence])
            
            if all_combinations:
                heatmap_df = pd.concat([heatmap_df, pd.DataFrame(all_combinations, columns=['Substrate', 'Enzyme', 'Percentage', 'Sequence'])], ignore_index=True)

            # Pivot table with explicit enzyme order
            heatmap_pivot = heatmap_df.pivot(index='Enzyme', columns='Substrate', values='Percentage')
            heatmap_pivot = heatmap_pivot.reindex(enzyme_list)
            
            # Prepare hover text matrix
            hover_texts = heatmap_df.apply(lambda row: f"{row['Enzyme']} - {row['Substrate']} ({row['Sequence']}): {row['Percentage']:.2f}%", axis=1)
            hover_text_matrix = heatmap_pivot.copy()
            for i, row in heatmap_df.iterrows():
                hover_text_matrix.at[row['Enzyme'], row['Substrate']] = hover_texts.iloc[i]

            fig_heatmap = px.imshow(
                heatmap_pivot, 
                labels=dict(x="Substrate", y="Enzyme", color="Percentage (%)"),
                x=substrate_list, 
                y=enzyme_list,
                aspect="auto", 
                color_continuous_scale=st.session_state.color_scale,
                zmin=0, 
                zmax=100
            )

            fig_heatmap.update_traces(
                hovertemplate="<b>Enzyme:</b> %{y}<br><b>Substrate:</b> %{x}<br><b>Sequence:</b> %{customdata}<br><b>Percentage:</b> %{z:.2f}%",
                customdata=hover_text_matrix.values
            )
            fig_heatmap.update_layout(
                coloraxis_colorbar=dict(
                    title="Percentage (%)",
                    tickvals=[0, 50, 100],
                    ticktext=["0%", "50%", "100%"]
                )
            )
            st.subheader("Heatmap of Total Light Emission (Percentage)")
            st.plotly_chart(fig_heatmap, use_container_width=True, key="heatmap")

# Function to display substrate structure images (unchanged)
def display_substrate_image(substrate_name):
    try:
        substrate_number = int(substrate_name.split()[-1])
        image_path = f"S{substrate_number}image.jpg"
        st.image(image_path, caption=f"Structure of {substrate_name}", use_container_width=True)
    except ValueError:
        st.error(f"Could not determine substrate number from '{substrate_name}'")

# Main logic to switch between welcome window and app
if not st.session_state.show_main_app:
    show_welcome_window()
else:
    show_main_app()

# Dropdown to select a substrate for image display (commented out)
# if st.session_state.show_main_app:
#     selected_substrate = st.selectbox("Select a substrate to view its structure:", substrate_list, index=0)
#     if selected_substrate:
#         display_substrate_image(selected_substrate)