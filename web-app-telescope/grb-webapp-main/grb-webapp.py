##############################################################################
# Imports
##############################################################################



import streamlit as st
import os
import shutil
import glob
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objs as go
from streamlit_plotly_events import plotly_events
from json import loads  
from astropy.coordinates import SkyCoord

import grblc.fitting.io as io
from grblc.fitting import Model, Lightcurve
import gcn_altered.scraper as scraper
import ads
import ads_altered
import re
from zipfile import ZipFile

# ignores UserWarning for non-GUI plotting
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")



##############################################################################
# Page config
##############################################################################



apptitle = "GRB Optical Repository"

st.set_page_config(page_title=apptitle, layout='centered')

st.title("Gamma Ray Bursts Optical Repository")



##############################################################################
# App contents
##############################################################################



pages = [
    "GRB Quickview",
    "Fit",
    "Bulk Download",
    "Information",
    "App Guide"
]
with st.sidebar.container():
    page_nav = st.selectbox("Menu", pages)

    
## GRB Quickview
##############################################################################


### Cleaning download directory on each new session

if 'reject_GRB' not in st.session_state:
    st.session_state.reject_GRB = 0

if 'event_GRB' not in st.session_state:
    st.session_state.event_GRB = 0

if 'count' not in st.session_state:            
    st.session_state.count = 0
        
if st.session_state.count == 0:
    if os.path.exists('gcncc/output/') == True:
        shutil.rmtree('gcncc/output')

if 'filenameBULK' not in st.session_state: 
    if os.path.exists('BD/') == True:
        shutil.rmtree('BD/')

### Loading data

@st.experimental_memo  # atlernative to @st.cache for dataframes
def load_data():
    ''' Reads the list of GRBs'''
    data =  pd.read_csv("mainsheet.csv", sep=",", engine="python", header=0, index_col=0, encoding="ISO-8859-1")
    return data
data = load_data()
grbexpander = st.sidebar.expander("GRB optical Repository",expanded=True)
select_event = grbexpander.selectbox("Find GRB by name", data.index)

### Raw data

if page_nav == "GRB Quickview":

    st.header("GRB " + str(select_event) + " Quickview")

    ra,dec,z,beta = st.columns(4)
    #z,beta = st.columns(2)
    ra.metric("Right Ascension",data.loc[data.index == select_event, "RA"].to_numpy()[0])
    dec.metric("Declination",data.loc[data.index == select_event, "DEC"].to_numpy()[0])
    z.metric("Redshfit",data.loc[data.index == select_event, "z"].to_numpy()[0])
    beta.metric("Spectral index", str(data.loc[data.index == select_event, "beta"].to_numpy()[0]) + ' ' + u"\u00B1" + ' ' +str(data.loc[data.index == select_event, "beta_err"].to_numpy()[0]))
    

    try:
        filename = data.loc[data.index == select_event, "path_F"].to_numpy()[0]  
               # matching the event name to the relavant file name
        #light = io.read_data(filename, header=None)
        light = pd.read_csv(filename, sep=r"\t+|\s+", engine="python", header=None, index_col=0)
        fig = px.scatter(
            data_frame=light,
            x=np.log10(light.index.values),
            y=np.log10(light[1].values),
            error_y= light[2].values / (light[1].values * np.log(10)),
            color=light[3].values,
            hover_data=[4],
            #title="Raw Data",
        )
        fig.update_layout(
            xaxis_title='log Time (s)',
            yaxis_title='log Flux (erg'+r'$cm^{-2} s^{-1}$'+')',
            font=dict(
                size=18
            )
        )
        st.plotly_chart(fig)
        #st.write("Source:", str(data.loc[data.index == select_event, "path"].values()[0]))
    except:
        filename = data.loc[data.index == select_event, "path_F"].to_numpy()[0]  
               # matching the event name to the relavant file name
        #light = io.read_data(filename, header=None)
        light = pd.read_csv(filename, sep=r"\t+|\s+", engine="python", header=0, index_col=0)
        fig = px.scatter(
            data_frame=light,
            x=np.log10(light.index),
            y=np.log10(light['flux'].values),
            error_y= light['flux_err'].values / (light['flux'].values * np.log(10)),
            color=light['band'].values,
            hover_data=['source'],
            #title="Raw Data",
        )
        fig.update_layout(
            xaxis_title='log Time (s)',
            yaxis_title='log Flux (erg'+r'$cm^{-2} s^{-1}$'+')',
            font=dict(
                size=18
            )
        )
        st.plotly_chart(fig)
        #st.write("Source:", str(data.loc[data.index == select_event, "path"].values()[0]))

    #### Download data
    
    mag, flux = st.columns(2)
    
    try:
        df = pd.read_csv(data.loc[data.index == select_event, "path_M"].to_numpy()[0], sep='\t', index_col=0)
        csv = df.to_csv(sep='\t').encode('utf-8')
        mag.download_button(
            label="Download data in magnitude",
            data= csv,
            file_name = data.loc[data.index == select_event, "path_M"].to_numpy()[0],
            mime='text/csv',
        )
    except (AttributeError,ValueError) as e:
        mag.warning("No magnitude file available")

    try:
        df = pd.read_csv(data.loc[data.index == select_event, "path_F"].to_numpy()[0], sep='\t', index_col=0)
        csv = df.to_csv(sep='\t').encode('utf-8')

        flux.download_button(
            label="Download data in erg "+r"$cm^{-2} s^{-1}$",
            data=csv,
            file_name = data.loc[data.index == select_event, "path_F"].to_numpy()[0],
            mime='text/csv',
        )
    except (AttributeError,ValueError) as e:
        flux.warning("No flux file available")


## Fitting
##############################################################################


elif page_nav == "Fit":

    
    st.header("Fit")

    ### Instructions

    instructions = st.expander("Instructions")
    instructions.write(r"""How to fit your selected data:
- Select the model for your GRB. For description of the models, read the "Information" page.
- Choose parameters and click 'Done'
- Next step is optional. You can reject data points by mouse clicks.
- Now, click 'Start' to see the 'Fit'.
- By default, all points are considered for fitting.
""")

    ### Choose model
    
    param = st.expander("Fit Model and Method")

    select_model = param.selectbox("Select model",
    ["None",
    "Willingale 2007",
    "Simple broken power law",
    "Smooth broken power law",
    ])

    with param.form("Select parameters"):

        if select_model == "Willingale 2007":

            model = Model.W07(vary_t=False)

            T = st.number_input(label="log Time", 
                max_value=10., 
                min_value=1e-10, 
                value=data.loc[data.index == select_event, 'logTa'].to_numpy()[0],
                )
            F = st.number_input(label="log Flux",  
                max_value=0., 
                min_value=-100., 
                value=data.loc[data.index == select_event, 'logFa'].to_numpy()[0],
                )
            alpha = st.number_input(label="alpha", 
                max_value=25.,
                min_value=-5.,
                value=data.loc[data.index == select_event, 'alpha'].to_numpy()[0],
                )
            t = 0.
            p = [T, F, alpha, t]

        elif select_model == "Simple broken power law":

            model = Model.SIMPLE_BPL()

            T = st.number_input(label="log Time", 
                max_value=10., 
                min_value=1e-10, 
                value=data.loc[data.index == select_event, 'T_bpl'].to_numpy()[0],
                )
            F = st.number_input(label="log Flux", 
                max_value=0., 
                min_value=-100., 
                value=data.loc[data.index == select_event, 'F_bpl'].to_numpy()[0],
                )
            alpha1 = st.number_input(label="alpha", 
                max_value=25.,
                min_value=-5.,
                value=data.loc[data.index == select_event, 'alpha1_bpl'].to_numpy()[0],
                )
            alpha2 = st.number_input(label="beta", 
                max_value=25.,
                min_value=0.,
                value=data.loc[data.index == select_event, 'alpha2_bpl'].to_numpy()[0],
                )
            p = [T, F, alpha1, alpha2]

        elif select_model == "Smooth broken power law":

            model = Model.SMOOTH_BPL()

            T = st.number_input(label="log Time", 
                max_value=10., 
                min_value=1e-10, 
                value=data.loc[data.index == select_event, 'T_smooth'].to_numpy()[0],
                )  
            F = st.number_input(label="log Flux",
                max_value=0., 
                min_value=-100., 
                value=data.loc[data.index == select_event, 'F_smooth'].to_numpy()[0],
                )
            alpha1 = st.number_input(label="alpha", 
                max_value=25.,
                min_value=-5.,
                value=data.loc[data.index == select_event, 'alpha1_smooth'].to_numpy()[0],
                )
            alpha2 = st.number_input(label="beta", 
                max_value=25.,
                min_value=0.,
                value=data.loc[data.index == select_event, 'alpha2_smooth'].to_numpy()[0],
                )
            S = st.number_input(label="S", 
                max_value=1.797e+308, # streamlit's inf
                min_value=-1.797e+308,
                value=data.loc[data.index == select_event, 'S'].to_numpy()[0],
                )
            p = [T, F, alpha1, alpha2, S]
            
        else: 
            model = Model.SIMPLE_BPL() # dummy model
            p = None

        select_mc = st.checkbox("Run MCMC")

        submit_model = st.form_submit_button("Save")

    ### Reject Points

    def sigma_reject(light, model, p):
        '''
        INPUT
        -----
        Data, Model, Parameters
        OUTPUT
        ------
        3 sigma rejected data
        '''

        flux_mod = model(np.array(light['time_sec']), *p)
        n = len(light)- len(p)
        light_temp = light.assign(res = np.array(light['flux'] - flux_mod))
        s = np.sqrt(np.sum([err ** 2 for err in light_temp['res'] if err is not None])/n)
        light_new = light_temp.loc[(light_temp.res < 3*s) & (light_temp.res > -3*s)]

        return light_new

    ### Make the fit

    def lcproduce(light, model):
        '''
        INPUT
        -----
        Data, Model, Parameters
        OUTPUT
        ------
        Light curve
        '''
        lc = Lightcurve(
            xdata=light["time_sec"],
            ydata=light["flux"],
            yerr=light["flux_err"],
            model=model,
            data_space='lin',
            name = "GRB " + select_event
        )
        
        return lc

    try:
        filename = data.loc[data.index == select_event, "path_F"].to_numpy()[0]  
               # matching the event name to the relavant file name
        light = io.read_data(filename)

        with st.expander("Reject Points"):

            fig = go.FigureWidget([go.Scatter(
                    x=np.log10(light["time_sec"]),
                    y=np.log10(light["flux"]),
                    error_y=dict(
                        type='data',
                        array=light["flux_err"] / (light["flux"] * np.log(10)),
                        visible=True),
                    mode='markers')])
            fig.update_layout(
                xaxis_title = 'log Time (s)',
                yaxis_title = 'log Flux (erg'+r'$cm^{-2} s^{-1}$'+')'
            )
            scatter = fig.data[0]
            colors = ['#000000'] * len(light)
            scatter.marker.color = colors

            if 'reject_index' not in st.session_state:
                st.session_state.reject_index = []

            if st.session_state.reject_GRB != st.session_state.event_GRB:
                st.session_state.reject_index = []
                st.session_state.reject_GRB = st.session_state.event_GRB

            clicked_point = plotly_events(fig, key='light')

            if len(clicked_point) != 0:
                st.session_state.reject_index.append(clicked_point[0]['pointIndex'])

            table_reject = go.FigureWidget([go.Table(
                header=dict(values=['time_sec','flux','flux_err','band'],
                            fill = dict(color='#eb6234'),
                            align = ['left'] * 5),
                cells=dict(values=[light[col][st.session_state.reject_index] 
                                for col in ['time_sec','flux','flux_err','band']],
                        fill = dict(color='#eff2f5'),
                        align = ['left'] * 5))])
            fig.update_layout(margin=dict(l=0,r=0,b=0,t=0))

            st.write("REJECTED POINTS:")
            st.plotly_chart(table_reject)

        if st.button("Start"):

            light = light.drop(index=st.session_state.reject_index, axis=0)

            if select_model == "None":
                st.warning("No model was selected")
                lc = lcproduce(light, model)
                st.set_option('deprecation.showPyplotGlobalUse', False)
                st.pyplot(lc.show_data(fig_kwargs=dict(dpi=72)))
            else:  
                try:
                    with st.spinner("Loading..."):
                        lc = lcproduce(light, model)
                        #light = sigma_reject(light, model,p) #will be added after initial guess added

                        if st.button('Stop',key='stop'):
                            raise KeyboardInterrupt

                        fit = st.empty()
                        chisq = st.empty()
                        corner = st.empty()
                        details = st.empty()

                        tt = data.loc[data.index == select_event, 'tt'].to_numpy()[0]
                        tf = data.loc[data.index == select_event, 'tf'].to_numpy()[0]
                        lc.set_bounds(xmin=tt, xmax=tf)

                        details.write(lc.fit(p0=p, run_mcmc=select_mc))
                        
                        st.set_option('deprecation.showPyplotGlobalUse', False)
                        fit.pyplot(lc.show_fit(detailed=False, print_res=False))
                        chisq.pyplot(lc.show_fit(show_plot=False, show_chisq=True, print_res=False))
                        if select_mc:
                            corner.pyplot(lc.show_fit(show_plot=False, show_chisq=False, show_corner=True, print_res=False))

                except KeyboardInterrupt:
                    st.stop()

        else:
            st.warning("Submit a model!")
        if st.button("Restart"):
            st.session_state.reject_index=[]
            st.legacy_caching.clear_cache()
            st.experimental_rerun()
    except:
        st.error('Fit needs flux converted files')



## GCN and ADS Search
#################################################################################################
elif page_nav == "GCN and ADS Search":

    st.header("Gamma-ray Coordinates Network(GCN) and NASA ADS Search")
    GRB_ID = st.text_input("Enter GRB ID",'-')

    #contain = st.container()
    col1, col2= st.columns(2)
    if GRB_ID != '-':
        col1.subheader('GCN results')
        col2.subheader('NASA ADS results')
    try : 
        if GRB_ID != '-':
            if not os.path.exists('gcncc/output/'+GRB_ID+'/'+GRB_ID+'_sentences_mag.txt'):
                scraperOBJ=scraper.Scraper()
                print('Collecting GCN circulars for GRB'+GRB_ID+'...')
                scraperOBJ.index_live_lookup(GRB_ID,streamlit_PB = True,container = col1)
                st.session_state.count += 1
            col1.dataframe(pd.read_csv('gcncc/output/'+GRB_ID+'/'+GRB_ID+
                           '_sentences_mag.txt',delimiter = '\t'),width = 400)
            print('Collected GCN circulars!...File ready to download')   
            with open('gcncc/output/'+GRB_ID+'/'+GRB_ID+'_sentences_mag.txt', 'rb') as f:
                if col1.download_button(label="Download GCN magnitude table", 
                                        data = f, 
                                        file_name=GRB_ID+'_sentences_mag.txt'):
                    col1.write('File Downloaded!')
    except (FileNotFoundError,ZeroDivisionError,TypeError) as e:
        col1.error("There is no GCN data available for the object '"+GRB_ID+
                   "' found! \nPlease confirm the input id")   
    
    
    
    try:
        if (GRB_ID != '-') & (len(GRB_ID)>=6):
            ads.config.token = 'AH0q6nzE3RfgbJRSGp7hf5X3CbLvriPqTcVCFhtU'
            dictonary = ads_altered.search.litSearch(GRB_ID,keywords = 'GRB')

            ADSdata = pd.DataFrame(columns=['Title','ADS link'])
            gap = 1/len(dictonary['articlelist'])
            progbar = col2.progress(0)
            for i,paper in enumerate(dictonary['articlelist']):
                ADSdata = ADSdata.append({'Title':paper.title[0], 
                                          'ADS link': 'https://ui.adsabs.harvard.edu/abs/'+paper.bibcode+'/abstract'},
                                          ignore_index=True)
                progbar.progress(gap+gap*i)

            col2.dataframe(ADSdata, width = 400)
            @st.cache
            def convert_df(df):
             # IMPORTANT: Cache the conversion to prevent computation on every rerun
             return df.to_csv().encode('utf-8')

            ADScsv = convert_df(ADSdata)

            col2.download_button(
                 label="Download table as CSV",
                 data=ADScsv,
                 file_name='ADSlinks.csv',
                 mime='text/csv',
             )
            
        elif (GRB_ID != '-'):
            col2.error("There is no NASA ADS data available for the object '"+GRB_ID+
                       "' found! \nPlease confirm the input id")
    except (ZeroDivisionError,TypeError,ConnectionError) as e:
        if e is ConnectionError:
            st.error("Connection Error! Try refreshing the session.")
        else : 
            col2.error("There is no NASA ADS data available for the object '"+GRB_ID+
                       "' found! \nPlease confirm the input id")
            
            
            
            
## Bulk Download
#############################################################################


elif page_nav == "Bulk Download":            

    def relationCHECKER(RAmax, RAmin, DECmax, DECmin, Zmax, Zmin):
        if (RAmax>RAmin) & (DECmax>DECmin) & (Zmax>Zmin):
            return 1
        else:
            return 0

    def download_option(RAmax, RAmin, DECmax, DECmin, Zmax, Zmin):
            if (np.sum(np.array([RAmax,RAmin,DECmax,DECmin,Zmax,Zmin]) == '-') >0) and not (relationCHECKER(RAmax,RAmin,DECmax,DECmin,Zmax,Zmin)):
                st.write("Enter proper values!")
            else:
                downloadDATA = data.loc[(data.RA_deg>RAmin) & (data.RA_deg<RAmax) & (data.DEC_deg>DECmin) & (data.DEC_deg<DECmax) & (data.z>=Zmin) & (data.z<=Zmax)]
                LCfiles = glob.glob('mastersample/*.txt')
                filenameBULK = 'BD/'+str(np.round(RAmin,decimals=1)) + str(np.round(RAmax,decimals=1)) + '_'  + str(np.round(DECmin,decimals=1)) + str(np.round(DECmax,decimals=1)) + '_' + str(np.round(Zmin,decimals=1)) + str(np.round(Zmax,decimals=1)) + '.zip'
                
                try:
                    os.mkdir('BD')
                except FileExistsError as e:
                    pass
                
                with ZipFile(filenameBULK, 'w') as zipObj2:
                    for i,row in downloadDATA.iterrows():
                        zipObj2.write(row.path,)

                return filenameBULK, len(downloadDATA)

    c = SkyCoord(data.RA.to_numpy()[0],data.DEC.to_numpy()[0],frame='icrs')

    data['RA_deg'] = c.ra.deg
    data['DEC_deg'] = c.dec.deg
    min_z = min(data.z.to_numpy())
    max_z = max(data.z.to_numpy())
    st.header('Bulk download LC data')
    
    if 'filename' not in st.session_state:
        st.session_state.filenameBULK = None
    
    with st.form("my_form"):
        RAmin,RAmax = st.slider(
             'Right Ascension',
             0.0, 360.0, step=0.0001, value =(0.0, 360.0))

        DECmin,DECmax = st.slider(
             'Declination',
             -90.0, 90.0, step=0.0001, value =(-90.0, 90.0))
        Zmin,Zmax = st.slider(
             'Redshift',
             min_z, max_z, step=0.0001, value =(float(min_z), float(max_z)))
        submitted = st.form_submit_button("Generate Download Link")
        if submitted:
            st.session_state.filenameBULK, length = download_option(RAmax,RAmin,DECmax,DECmin,Zmax,Zmin)        
    
    
    if st.session_state.filenameBULK != None:
        with open(st.session_state.filenameBULK, 'rb') as f:
            st.write(str(length) + ' LC files found, Click the below button to download the zip file.')
            #st.write(st.session_state.filenameBULK)
            st.download_button('Download Zip', f, file_name=st.session_state.filenameBULK)
    
            
## Information
##############################################################################


elif page_nav == "Information":
    st.header("Information")

    st.write("### Simple broken law")
    st.write("""
This is an empirical piece-wise model for GRB lightcurve afterglows.
The function is as follows:
        """)
    st.latex(r'''
\begin{split}f(t) = 
\left \{ \begin{array}{ll} \displaystyle{F_i \left (\frac{t}{T_i} \right)^{-\alpha} } & {\rm for} \ \ t < T_i \\ 
\displaystyle{F_i \left ( \frac{t}{T_i} \right )^{-\beta} } & {\rm for} \ \ t \ge T_i, \\ \end{array} \right . \end{split}
        ''')
    st.write(r'''
where the transition from the exponential to the power law occurs at the point $T_i$, $F_i$,
$\alpha$ determines the temporal decay index of the initial power law, and 
$\beta$ is the temporal decay index of the final power law.
As implemented, log space is used for the time sec and flux erg cm$^{-2}$s$^{-1}$. 
This means that for a light curve in which the afterglow plateau phase ends at 10,000 seconds corresponds to a $T_i$ of 5.
        ''')


    st.write("### Smooth broken law")
    st.write("""
This is an empirical piece-wise model for GRB lightcurve afterglows.
The function is as follows:
        """)
    st.latex(r'''
f(t) = F_i \left (\left (\frac{t}{T_i} \right )^{S\alpha} + 
                \left (\frac{t}{T_i} \right )^{S \beta} \right )^{-\frac{1}{S}}
        ''')
    st.write(r'''
where the transition from the exponential to the power law occurs at the point $T_i$, $F_i$,
$\alpha$ determines the temporal decay index of the initial power law,  
$\beta$ is the temporal decay index of the final power law, and
$S$ is the smoothing factor
As implemented, log space is used for the time sec and flux erg cm$^{-2}$s$^{-1}$. 
This means that for a light curve in which the afterglow plateau phase ends at 10,000 seconds corresponds to a $T_i$ of 5.
        ''')


    st.write("### Willingale et al. (2007) model")
    st.write("""
    This is a phenomenological model for GRB lightcurve afterglows popularized in the paper by Willingale et. al, (2007).
    Taken from his paper, it is as follows:
        """)
    st.latex(r'''
\begin{split}f(t) = \left \{ \begin{array}{ll}\displaystyle{F_i \exp{\left ( \alpha_i 
                    \left( 1 - \frac{t}{T_i} \right) \right )} \exp{\left (- \frac{t_i}{t} \right )}} & {\rm for} \ \ t < T_i \\ ~ & ~ \\ 
                    \displaystyle{F_i \left ( \frac{t}{T_i} \right )^{-\alpha_i} \exp{\left ( - \frac{t_i}{t} \right )}} & {\rm for} \ \ t \ge T_i, \\\end{array} \right .\end{split}
        ''')
    st.write(r'''
where the transition from the exponential to the power law occurs at the point $T_i$, $F_i$,
$\alpha$ determines the temporal decay index of the initial power law, and
$t_i$ is the time of the initial rise of the lightcurve.
As implemented, log space is used for the time sec and flux erg cm$^{-2}$s$^{-1}$. 
This means that for a light curve in which the afterglow plateau phase ends at 10,000 seconds corresponds to a $T_i$ of 5.
        ''')


## App Guide
##############################################################################


elif page_nav == "App Guide":
    
    st.header("Here is a video demo")
    st.video('video.webm', format="video/webm")

    
    
    

    
