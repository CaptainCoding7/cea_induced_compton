from kivy.app import App
from kivy.uix.gridlayout import GridLayout
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.core.window import Window
from kompaneets_cc2 import mainKompaneetsFromGui

class showGUI(GridLayout):

    def __init__(self, **kwargs):
        super(showGUI, self).__init__(**kwargs)
        self.cols = 4
        self.lines = 5
        self.title = "test"
        ######################### Text fields of the physical parameters
        # kT
        self.add_widget(Label(text='kT (keV)'))
        self.kT = TextInput(text="80",multiline=False)
        self.add_widget(self.kT)
        # B
        self.add_widget(Label(text='B (G)'))
        self.B = TextInput(text="1e5",multiline=False)
        self.add_widget(self.B)
        #Thompson optical depth
        self.add_widget(Label(text='Thompson Optical Depth'))
        self.pT = TextInput(text="2.36",multiline=False)
        self.add_widget(self.pT)
        # Schwarzschild radius
        self.add_widget(Label(text='Schwarzschild radius (cm)'))
        self.R = TextInput(text="5e7",multiline=False)
        self.add_widget(self.R)
        
        ###################### Graph parameters
        self.add_widget(Label(text='xmin (keV)'))
        self.xmin = TextInput(text="1e-5",multiline=False)
        self.add_widget(self.xmin)
                
        self.add_widget(Label(text='xmax (keV)'))
        self.xmax = TextInput(text="1e4",multiline=False)
        self.add_widget(self.xmax)

        self.add_widget(Label(text='ymin (erg/s)'))
        self.ymin = TextInput(text="1e28",multiline=False)
        self.add_widget(self.ymin)
                
        self.add_widget(Label(text='ymax (erg/s)'))
        self.ymax = TextInput(text="5e33",multiline=False)
        self.add_widget(self.ymax)
        
        
        self.btn = Button(text="PLOT")
        self.btn.bind(on_press=self.buttonClicked)
        self.add_widget(self.btn)
        
        # button click function
    def buttonClicked(self,btn):
        
        mainKompaneetsFromGui(float(self.R.text),float(self.B.text),float(self.kT.text),float(self.pT.text),
                              float(self.xmin.text),float(self.xmax.text),float(self.ymin.text),float(self.ymax.text))
        
class MyApp(App):

    def build(self):
        self.title = "Plotting parameters"
        Window.size = (800, 150)
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