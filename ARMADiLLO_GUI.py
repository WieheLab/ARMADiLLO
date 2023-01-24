#!/usr/bin/python3

#libraries
import xlrd, re
import glob
import os
import string
import sys
import argparse

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk,Gdk

class GuiHandler:
    def __init__(self):
        self.singleFile=""
        self.seqFile=""
        self.ucaFile=""
        self.outFasta=""
        self.outputFormat="HTML"
        self.species="human"
        self.chain="heavy"
        
        builder=Gtk.Builder()
        builder.add_from_file("armadillo_GUI.glade")
        builder.connect_signals(self)
        
        self.window=builder.get_object("window")
        self.window.connect("destroy",Gtk.main_quit)

        self.file_selectors=[]
        self.file_selectors.append(builder.get_object("SingleFile"))
        self.file_selectors.append(builder.get_object("SeqFile"))
        self.file_selectors.append(builder.get_object("UCAFile"))
        
        self.fileLabel=builder.get_object("CountFile")
        self.fileLabel.set_label("N/A")

        css=b'''
        button {
        background:white;
        border-radius: 20px;
        }
        '''

        style_provider=Gtk.CssProvider()
        style_provider.load_from_data(css)
        Gtk.StyleContext.add_provider_for_screen(
            Gdk.Screen.get_default(),
            style_provider,
            Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION
            )

        
        self.window.show_all()


    def on_RunButton_clicked(self,widget):

        command="./ARMADiLLO -m Mutability.csv -s Substitution.csv -max_iter 10000 "
        if self.outputFormat=="HTML":
            command+="-HTML "
        elif self.outputFormat=="Simple Text":
            command+="-simple_text "
        elif self.outputFormat=="Standard Text":
            command+="-text "
        elif self.outputFormat=="Full Output":
            command+="-fulloutput "
        else:
            command+="-fulloutput "

        if self.species=="Human":
            command+="-species human "
        elif self.species=="Mouse":
            command+="-species mouse "
        elif self.species=="Rhesus":
            command+="-species rhesus "
        else:
            command+=""

        if self.chain=="Heavy":
            command+="-chain heavy "
        elif self.chain=="Lambda":
            command+="-chain lambda "
        elif self.chain=="Kappa":
            command+="-chain kappa "
        else:
            command+=""

        if self.singleFile:
            if "yaml" in self.singleFile or "cvs" in self.singleFile:
                command+="-partis {}".format(self.singleFile)
            elif "fasta" in self.singleFile:
                command+="-SMUA {}".format(self.singleFile)
            else:
                print("problem with selecting single file")
        elif self.seqFile and self.ucaFile:
            command+="-seq {} -uca {}".format(self.seqFile,self.ucaFile)
        else:
            print("no files selected")
            
        if self.fileLabel.get_label()!="N/A":
            command+=" > {}".format(self.fileLabel.get_label())

        os.system(command)
        print("need a done command")
        
    def on_ExitButton_clicked(self,widget):
        Gtk.main_quit()

    def on_SingleFile_file_set(self,widget):
        self.singleFile=widget.get_filename()
        if self.fileLabel.get_label()=="N/A":
            outfile=os.path.splitext(self.singleFile)[0]+".txt"
            self.fileLabel.set_label(os.path.basename(outfile))
    def on_SeqFile_file_set(self,widget):
        self.seqFile=widget.get_filename()
        if self.fileLabel.get_label()=="N/A":
            outfile=os.path.splitext(self.seqFile)[0]+".txt"
            self.fileLabel.set_label(os.path.basename(outfile))
    def on_UCAFile_file_set(self,widget):
        self.ucaFile=widget.get_filename()

    def on_OutputHTMLRadio_toggled(self,widget):
        self.outputFormat=widget.get_child().get_label()

    def on_HumanRadio_toggled(self,widget):
        self.species=widget.get_child().get_label()

    def on_HeavyRadio_toggled(self,widget):
        self.chain=widget.get_child().get_label()

    def on_ClearFiles_clicked(self,widget):
        for filechooser in self.file_selectors:
            filechooser.unselect_all()
        self.fileLabel.set_label("N/A")

    def on_savefile_clicked(self,widget):
        dialog=Gtk.FileChooserDialog(title="File to save data to...",parent=self.window,action=Gtk.FileChooserAction.SAVE,buttons=(Gtk.STOCK_CANCEL,Gtk.ResponseType.CANCEL,Gtk.STOCK_OPEN,Gtk.ResponseType.OK))
        cansave=False
        response=dialog.run()
        if response == Gtk.ResponseType.OK:
            print(dialog.get_filename())
            self.outFasta=dialog.get_filename()
            
            if os.path.exists(self.outFasta) == True:  # does file already exists?
                dialog2 = DialogSaveFile(self.window, self.outFasta)  # ask to confirm overwrite
                response2 = dialog2.run()
                if response2 == Gtk.ResponseType.OK:
                    cansave = True
                else:
                    self.outFasta=""
                dialog2.destroy()
            else:
                cansave = True
            self.fileLabel.set_label(self.outFasta)

        elif response == Gtk.ResponseType.CANCEL:
            print("Canceled selection")

        dialog.destroy()
        
        
class DialogSaveFile(Gtk.Dialog):
    def __init__(self, parent, db):
        Gtk.Dialog.__init__(self, "Confirm overwrite", parent, 0,
            (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
             Gtk.STOCK_OK, Gtk.ResponseType.OK))
        self.box = self.get_content_area()
        self.label = Gtk.Label("The file `" + db + "` exists.\nDo you want it to be overwritten?")
        self.box.add(self.label)
        self.show_all()

        
        
def main():
    g=GuiHandler()
    Gtk.main()
    return
    

if __name__=="__main__":
    main()
    

