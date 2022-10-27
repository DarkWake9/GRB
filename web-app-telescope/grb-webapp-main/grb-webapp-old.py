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
import plotly.offline as py
py.init_notebook_mode()
from streamlit_plotly_events import plotly_events
from json import loads  
from astropy.coordinates import SkyCoord

import grblc
import io_altered
import gcn_altered.scraper as scraper
import ads
import ads_altered
import re
from zipfile import ZipFile
##############################################################################
# Workaround for the limited multi-threading support in matplotlib.
# Per the docs, we will avoid using `matplotlib.pyplot` for figures:
# https://matplotlib.org/3.3.2/faq/howto_faq.html#how-to-use-matplotlib-in-a-web-application-server.
# Moreover, we will guard all operations on the figure instances by the
# class-level lock in the Agg backend.
##############################################################################
from matplotlib.backends.backend_agg import RendererAgg
_lock = RendererAgg.lock



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
    "GCN and ADS Search",
    "Bulk Download",
    "Information",
    "App Guide"
]
with st.sidebar.container():
    page_nav = st.selectbox("Menu", pages)



if 'count' not in st.session_state:            ### for cleaning GCN download directory on each new session
        st.session_state.count = 0
        
if st.session_state.count == 0:
    if os.path.exists('gcncc/output/') == True:
        shutil.rmtree('gcncc/output')

if 'filenameBULK' not in st.session_state:    ### for cleaning Bulk download directory on each new session
    if os.path.exists('BD/') == True:
        shutil.rmtree('BD/')
    
### Loading data

@st.experimental_memo  # atlernative to @st.cache for dataframes
def load_data():
    ''' Reads the list of GRBs'''
    data = pd.read_csv(
        "grb-list.csv", delimiter=","
    )  # list of 180 grb from master sample
    return data
data = load_data()
grbexpander = st.sidebar.expander("GRB optical Repository",expanded=True)
select_event = grbexpander.selectbox("Find GRB by name", data["GRB"])
filename = data.loc[data.GRB == select_event, "path"].to_numpy()[0]  
               # matching the event name to the relavant file name
light = io_altered.read_data(filename)





####################### Bulk Download functions
######################################################################


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


        
        


### Raw data
############################################################################

if page_nav == "GRB Quickview":

    st.header("GRB " + str(select_event) + " Quickview")

    ra,dec,z = st.columns(3)
    ra.metric("Right Ascension",data.loc[data.GRB == select_event, "RA"].to_numpy()[0])
    dec.metric("Declination",data.loc[data.GRB == select_event, "DEC"].to_numpy()[0])
    z.metric("Redshfit",data.loc[data.GRB == select_event, "z"].to_numpy()[0])
    
    fig = px.scatter(
        data_frame=light,
        x=np.log10(light["time_sec"]),
        y=np.log10(light["flux"]),
        error_y= light["flux_err"] / (light["flux"] * np.log(10)),
        color=light["band"],
        title="Raw Data",
    )
    fig.update_layout(
        xaxis_title='log Time (s)',
        yaxis_title='log Flux (erg cm^{-2} s^{-1})'
    )
    st.plotly_chart(fig)


    #### Download data
    @st.experimental_memo
    def convert_df(df):
        # IMPORTANT: Cache the conversion to prevent computation on every rerun
        return df.to_csv().encode('utf-8')
    csv = convert_df(light)
    st.download_button(
         label="Download data as CSV",
         data=csv,
         file_name=select_event+'_'+
                   data.loc[data.GRB == select_event, "datatype"].to_numpy()[0]+'.csv',
         mime='text/csv',
     )


## Fitting
##############################################################################


elif page_nav == "Fit":

    
    st.header("Fit")

    ### Instructions

    instructions = st.expander("Instructions")
    instructions.write(r"""How to obtain light curve fit:
- Select the model.
- Choose parameters and 'Submit'.
- Next step is optional. You can reject data points by mouse clicks.
- Once done rejecting, click 'Done' to see the 'Fit'.
- By default, all points are considered for fitting.
Refer 'Information' page to learn more about the models.
""")

    ### Choose model

    models = [
    "None",
    "Willingale 2007",
    "Simple broken power law",
    "Smooth broken power law",
    ]  # models available in the grblc package
    select_model = st.selectbox("Select model", models)
    para = st.expander("Select parameters")

    with para.form("Select parameters"):

        if select_model == "Willingale 2007":
            T = st.number_input(label="log Time", 
                min_value=1e-10, 
                max_value=10., 
                value=5.,
                )  # min, max, default # need to adjust
            F = st.number_input(label="log Flux",  
                min_value=-100., 
                max_value=0., 
                value=-12.5,
                )
            alpha = st.number_input(label="alpha", 
                min_value=-25.,
                max_value=5.,
                value=-0.1,
                )
            t = 0.
            p = [T, F, alpha, t]

        elif select_model == "Simple broken power law":
            T = st.number_input(label="log Time", 
                min_value=1e-10, 
                max_value=10., 
                value=5.,
                )  # min, max, default # need to adjust
            F = st.number_input(label="log Flux", 
                min_value=-100., 
                max_value=0., 
                value=-12.5,
                )
            alpha1 = st.number_input(label="alpha", 
                min_value=-25.,
                max_value=5.,
                value=-0.1,
                )
            alpha2 = st.number_input(label="beta", 
                min_value=-25.,
                max_value=0.,
                value=-1.5,
                )
            p = [T, F, alpha1, alpha2]

        elif select_model == "Smooth broken power law":
            T = st.number_input(label="log Time", 
                min_value=1e-10, 
                max_value=10., 
                value=5.,
                )  # min, max, default # need to adjust
            F = st.number_input(label="log Flux",
                min_value=-100., 
                max_value=0., 
                value=-12.5,
                )
            alpha1 = st.number_input(label="alpha", 
                min_value=-25.,
                max_value=5.,
                value=-0.1,
                )
            alpha2 = st.number_input(label="beta", 
                min_value=-25.,
                max_value=0.,
                value=-1.5,
                )
            S = st.number_input(label="S", 
                min_value=-100.,
                max_value=100.,
                value=1.,
                )
            p = [T, F, alpha1, alpha2, S]
            
        else: 
            pass

        select_mc = st.checkbox("Run MCMC")
        submit_model = st.form_submit_button("Submit")

    @st.experimental_memo
    def lcproduce(light, select_model, p):
        '''
        INPUT
        -----
        Data, Model, Parameters
        OUTPUT
        ------
        Light curve
        '''

        if select_model == "Willingale 2007":
            model = grblc.Model.W07(vary_t=False)

            lc = grblc.Lightcurve(
                xdata=light["time_sec"],
                ydata=light["flux"],
                yerr=light["flux_err"],
                model=model,
                data_space='lin',
            )
            return lc

        elif select_model == "Simple broken power law":
            model = grblc.Model.SIMPLE_BPL()

            lc = grblc.Lightcurve(
                xdata=light["time_sec"],
                ydata=light["flux"],
                yerr=light["flux_err"],
                model=model,
                data_space='lin',
            )
            return lc

        elif select_model == "Smooth broken power law":
            model = grblc.Model.SMOOTH_BPL()

            lc = grblc.Lightcurve(
                xdata=light["time_sec"],
                ydata=light["flux"],
                yerr=light["flux_err"],
                model=model,
                data_space='lin',
            )
            return lc

        else: 
            return None

    ### Reject Points

    with st.expander("Reject Points"):
        fig = go.FigureWidget([go.Scatter(
                x=np.log10(light["time_sec"]),
                y=np.log10(light["flux"]),
                error_y=dict(
                    type='data', # value of error bar given in data coordinates
                    array=light["flux_err"] / (light["flux"] * np.log10(10)),
                    visible=True),
                mode='markers')])
        fig.update_layout(
            xaxis_title = 'log Time (s)',
            yaxis_title = 'log Flux (erg cm^{-2} s^{-1})'
        )
        scatter = fig.data[0]
        colors = ['#000000'] * len(light)
        scatter.marker.color = colors

        if 'reject_index' not in st.session_state:
            st.session_state.reject_index = []

        clicked_point = plotly_events(fig, key='light')

        if len(clicked_point) != 0:
            i = clicked_point[0]['pointIndex']
            st.session_state.reject_index.append(i)

        table_reject = go.FigureWidget([go.Table(
            header=dict(values=['time_sec','flux','flux_err','band'],
                        fill = dict(color='#eb6234'),
                        align = ['left'] * 5),
            cells=dict(values=[light[col][st.session_state.reject_index] 
                               for col in ['time_sec','flux','flux_err','band']],
                       fill = dict(color='#eff2f5'),
                       align = ['left'] * 5))])
        fig.update_layout(margin=dict(l=0,r=0,b=0,t=0))
        
        st.plotly_chart(table_reject)

    ### Make the fit

    if st.button('Done', key='done'):
        light = light.drop(index=st.session_state.reject_index, axis=0)
        if select_model == "None":
            st.error("""No fit for 'None'""")
        else:  
            try:
                with st.spinner("Loading..."):
                    lc = lcproduce(light, select_model, p)
                    if st.button('Stop',key='stop'):
                        raise KeyboardInterrupt
                    with _lock:
                        fit = st.empty()
                        chisq = st.empty()
                        corner = st.empty()
                        st.write(lc.fit(p0=p, run_mcmc=select_mc))
                        st.set_option('deprecation.showPyplotGlobalUse', False)
                        fit.pyplot(lc.show_fit(detailed=False, print_res=False))
                        chisq.pyplot(lc.show_fit(show_chisq=True, print_res=False))
                        if select_mc:
                            corner.pyplot(lc.show_fit(show_corner=True, print_res=False))
                        
            except KeyboardInterrupt:
                st.stop()
        del st.session_state.reject_index
    else:
        st.warning("Submit a model!")



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
                scraperOBJ.grb_live_lookup(GRB_ID,streamlit_PB = True,container = col1)
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
    c = SkyCoord(data.RA.to_numpy(),data.DEC.to_numpy(),frame='icrs')

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


## Fitting
##############################################################################


elif page_nav == "App Guide":
    
    st.header("Here is a video demo")
    st.video('video.webm', format="video/webm")

    
    
    

    
