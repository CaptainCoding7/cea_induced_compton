from kivy.app import App
from kivy.uix.gridlayout import GridLayout
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.uix.checkbox import CheckBox
from kivy.core.window import Window
from kompaneets_cc2 import mainKompaneetsFromGui
from numpy import empty

class showGUI(GridLayout):

    def __init__(self, **kwargs):
        super(showGUI, self).__init__(**kwargs)
        self.cols = 4
        self.lines = 6
        self.title = "test"
        
        modelFile = open("param_cygX-1.txt")
        #modelFile = open("param_gx339-4.txt")
        mdlParam = modelFile.readlines()

        ######################### Text fields of the physical parameters
        # kT
        self.add_widget(Label(text='kT (keV)'))
        self.kT = TextInput(text=mdlParam[1].strip(),multiline=False)
        self.add_widget(self.kT)
        # B
        self.add_widget(Label(text='B (G)'))
        self.B = TextInput(text=mdlParam[2].strip(),multiline=False)
        self.add_widget(self.B)
        #Thompson optical depth
        self.add_widget(Label(text='Thomson Optical Depth'))
        self.pT = TextInput(text=mdlParam[3].strip(),multiline=False)
        self.add_widget(self.pT)
        # Corona radius
        self.add_widget(Label(text='Corona radius (cm)'))
        #self.label.bind(size=self.label.setter('text_size')    
        self.R = TextInput(text=mdlParam[4].strip(),multiline=False)
        self.add_widget(self.R)
        
        ###################### Graph parameters
        self.add_widget(Label(text='xmin (keV)'))
        self.xmin = TextInput(text="1e-5",multiline=False)
        self.add_widget(self.xmin)
                
        self.add_widget(Label(text='xmax (keV)'))
        self.xmax = TextInput(text="1e4",multiline=False)
        self.add_widget(self.xmax)

        self.add_widget(Label(text='ymin (erg/s)'))
        self.ymin = TextInput(text="1e22",multiline=False)
        self.add_widget(self.ymin)
                
        self.add_widget(Label(text='ymax (erg/s)'))
        self.ymax = TextInput(text="5e38",multiline=False)
        self.add_widget(self.ymax)
        
        self.add_widget(Label(text='Plot synchrotron/Brems contribution'))
        self.checkbox1 = CheckBox(active=True)
        self.plotSync = True
        self.add_widget(self.checkbox1)
        self.checkbox1.bind(active=self.on_checkbox1_active)

        self.add_widget(Label(text='Plot Compton contribution'))
        self.checkbox2 = CheckBox(active=True)
        self.plotComp = True
        self.add_widget(self.checkbox2)
        self.checkbox2.bind(active=self.on_checkbox2_active)
        
        self.add_widget(Label(text='Plot without induced Compton'))
        self.checkbox3 = CheckBox(active=True)
        self.plotIC = True 
        self.add_widget(self.checkbox3)
        self.checkbox3.bind(active=self.on_checkbox3_active)
        
        
        self.btn = Button(text="PLOT")
        self.btn.bind(on_press=self.buttonClicked)
        self.add_widget(self.btn)
        
        # button click function
    def buttonClicked(self,btn):
        
        mainKompaneetsFromGui(float(self.R.text),float(self.B.text),float(self.kT.text),float(self.pT.text),
                              float(self.xmin.text),float(self.xmax.text),float(self.ymin.text),float(self.ymax.text),
                              self.plotIC,self.plotSync,self.plotComp)

    
    def on_checkbox1_active(self,checkbox, value):
        self.plotSync = value
    
    def on_checkbox2_active(self,checkbox, value):
        self.plotComp = value
         
    def on_checkbox3_active(self,checkbox, value):
        self.plotIC = value   
        
            
class MyApp(App):

    def build(self):
        self.title = "Plotting parameters"
        Window.size = (1000, 200)
        Window.bind(on_request_close=self.close_app)
        return showGUI()

    # method which will render our application
    def close_app(self, *argss):
        # closing application
        App.get_running_app().stop()
        # removing window
        Window.close()
        
        
if __name__ == '__main__':
    MyApp().run()