#!/usr/bin/env python

from Tkinter import *
import ttk
import tkFileDialog
import os
import string
import os.path
import subprocess
import skmfunc
import sys
import shutil
import random

mycwd=string.join(sys.argv[0].split("\\")[:-1],"/")
ctlfile=open(mycwd+"/PhyTkPi.ctl").readlines()
for i in ctlfile:
    if "mafft_loc" in i:
        mafft_loc=i.split("=")[1][:-1]
    if "RAxML_loc" in i:
        RAxML_loc=i.split("=")[1][:-1]
    if "FigTree_loc" in i:
        FigTree_loc=i.split("=")[1][:-1]


def filedialog(*args):
    fili= tkFileDialog.askopenfilename(parent=root, title='choose a file')
    if fili != None:
        fastafile=open(fili).read()
        fastadisplay.delete(1.0,'end')
        fastadisplay.insert(1.0,fastafile)
        file_loc.set('File location: '+fili)
        seqnum.set('number of sequences: '+str(string.count(fastafile,">")))
        dt='Nucleotide'
        for i in ['q','e','i','o','p','f','j','l']:
            for k in string.split(fastafile,'\n'):
                if '>' in k:
                    pass
                else:
                    if i in string.lower(k):
                        dt="Protein"
        datatype.set('Data type: '+dt)
        for i in [il1,il2,il3]:
            i.configure(borderwidth=2,relief='sunken')
            i.grid_configure(padx=10,pady=10)
def updatealignopt(*args):
    if str(alignalg.get())=='MAFFT':
        mafftoptframer1.grid()
        mafftoptframer2.grid()
    elif str(alignalg.get())=='MUSCLE':
        mafftoptframer1.grid_remove()
        mafftoptframer2.grid_remove()
def mafftoptions(*args):
    if MAFFTalgorithm.get()=='accurate (1000 iterations)':
        pass
def alignmentgo(*args):
    fastaout=open(mycwd+'/tk.tempfiles/tempfasta.tempfasta','wb')
    fastaout.write(fastadisplay.get(1.0,'end'))
    fastaout.close()
    abatchout=open('C:/temp/tempbatchtk.bat','wb')
    abatchoutlist=[]
    if alignalg.get()=='MAFFT':
        abatchoutlist.append('"'+os.path.abspath(mycwd+'/tk.supportfiles/mafft.bat')+'"')
        if MAFFTalgorithm.get()=='Auto':
            pass
        elif MAFFTalgorithm.get()=='Fast (no iterations)':
            abatchoutlist.append(' --maxiterate 0')
    abatchoutlist.append(" '"+os.path.abspath(mycwd+'/tk.tempfiles/tempfasta.tempfasta')+'\' > "'+os.path.normpath(mycwd+'/tk.tempfiles/tempalign.temptempalign')+'"')
    abatchout.writelines(abatchoutlist)
    abatchout.close()
    subprocess.call('C:/temp/tempbatchtk.bat')
    os.remove('C:/temp/tempbatchtk.bat')
    tabalignment=string.split(skmfunc.fastatotab(mycwd+'/tk.tempfiles/tempalign.temptempalign'),'\n')
    for i in tabalignment[:]:
        tabalignment[tabalignment.index(i)]=string.split(i,'\t')
    tabalignment2=[]
    for i in tabalignment:
        if len(i)>1:
            tabalignment2.append(i)
    seqalign=[]
    namealign=[]
    for i in tabalignment2:
        seqalign.append(i[1])
        namealign.append(i[0])
    try: namealign.remove('\n')
    except ValueError: pass
    seqalignment.delete(1.0,'end')
    namealignment.delete(1.0,'end')
    seqalignment.insert(1.0,string.join(seqalign,'\n'))
    namealignment.insert(1.0,string.join(namealign,'\n'))
def runtreesearch(*args):
    print tbmethod.get()
    if tbmethod.get()=="RAxML":
        batch=open('C:/Temp/tempbatchtk.bat','w')
        raxmlopt=[RAxML_loc]
        raxmlopt.append(' -s "'+os.path.normpath(mycwd+'/tk.tempfiles/tempalign.temptempalign')+'" -m '+raxmlmodel.get()+' -w "'+os.path.normpath(mycwd+'/tk.tempfiles')+'" -f d -p '+str(int(random.uniform(0,999)))+' -n temptree.tre')
        batch.write(string.join(raxmlopt,''))
        batch.close()
        subprocess.call('C:/Temp/tempbatchtk.bat')
        batch=open('C:/Temp/tempbatchtk.bat','w')
        batch.write(FigTree_loc+" \""+os.path.normpath(mycwd+'/tk.tempfiles')+"/RAxML_bestTree.temptree.tre\"")
        batch.close()
        subprocess.call('C:/Temp/tempbatchtk.bat')
        os.remove(os.path.normpath(mycwd+'/tk.tempfiles')+"/RAxML_info.temptree.tre")
        #os.remove('C:/Temp/tempbatchtk.bat')


root = Tk()
root.title("Sean's phylogenetic pipeline")
root.columnconfigure(0,weight=1)

mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)

importframecolapse=True
if importframecolapse==True:
    importframe=ttk.Frame(mainframe, padding="3 3 12 12",width=1000)
    importframe.grid(column=0,row=1,sticky=(N,W,E,S))
    fviewframe=ttk.Frame(mainframe, padding="3 3 12 12", width=1000, height=300)
    fviewframe.grid(column=0,row=2,sticky=(N,W,E,S))
    alignoptframe=ttk.Frame(mainframe, padding="3 3 12 12",width=1000)
    alignoptframe.grid(column=0,row=3,sticky=(N,W,E,S))
    aviewframe=ttk.Frame(mainframe, padding="3 3 12 12",width=1000, height=300)
    aviewframe.grid(column=0,row=4,sticky=(N,W,E,S))
    treebuildframe=ttk.Frame(mainframe, padding="3 3 12 12",width=1000)
    treebuildframe.grid(column=0,row=5,sticky=(N,W,E,S))
    importframe.rowconfigure(0,weight=1)
    
    
    
    file_loc = StringVar()
    seqnum = StringVar()
    datatype= StringVar()
    
    icon=PhotoImage(file=mycwd+'\\bitmapico.gif')
    icondisplay=ttk.Label(importframe,image=icon,text='biorover productions',compound='bottom')
    icondisplay.grid(column=0,row=1,sticky=(N,E))
    importbutton=ttk.Button(importframe,text="Import fasta",command=filedialog)
    importbutton.grid(column=1,row=1,sticky=E)
    il1=ttk.Label(importframe,textvariable=file_loc)
    il1.grid(column=3,row=1,sticky=E)
    il2=ttk.Label(importframe,textvariable=seqnum)
    il2.grid(column=4,row=1,sticky=E)
    il3=ttk.Label(importframe,textvariable=datatype)
    il3.grid(column=5,row=1,sticky=E)
    
    fastadisplay=Text(fviewframe, height=10)    
    fdscrollh=ttk.Scrollbar(fviewframe,orient=HORIZONTAL,command=fastadisplay.xview)
    fdscrollh.grid(column=0,row=1,sticky=(N,S,E,W))
    fdscrollv=ttk.Scrollbar(fviewframe,orient=VERTICAL,command=fastadisplay.yview)
    fdscrollv.grid(column=1,row=0,sticky=(N,S,E,W))
    fastadisplay.configure(xscrollcommand=fdscrollh.set, yscrollcommand=fdscrollv.set,wrap='none')
    fastadisplay.grid(column=0,row=0,sticky=(N,S,E,W))

againiwanttocollapse=True
if againiwanttocollapse==True:
    alignalg=StringVar()
    alignalg.set('MAFFT')
    MAFFTalgorithm=StringVar()
    MAFFTalgorithm.set('Compromise (2 iterations)')
    MAFFToptions=StringVar()
    MAFFTdistance=StringVar()
    MAFFTdistance.set('6mers (faster)')
    MAFFTiterations=StringVar()
    MAFFTiterations.set('1000')
    acframe=ttk.Frame(alignoptframe)
    acframe.grid(row=1,column=1)
    alignchoose=ttk.Combobox(acframe,textvariable=alignalg,state='readonly',values=('MAFFT','MUSCLE'))
    alignchoose.grid(column=1,row=1,sticky=NE)
    alignchoose.bind('<<ComboboxSelected>>',updatealignopt)
    abframe=ttk.Frame(alignoptframe)
    abframe.grid(row=1,column=3,sticky=E)
    alignbutton=ttk.Button(abframe,text='Align',command=alignmentgo)
    alignbutton.grid(row=1,column=1)
    mafftoptframer1=ttk.Frame(alignoptframe)
    mafftoptframer1.grid(row=1,column=2)
    mafftoptframer2=ttk.Frame(alignoptframe)
    mafftoptframer2.grid(row=2,column=2,sticky=W)
    ttk.Label(mafftoptframer1,text='Alignment algorithm:').grid(row=1,column=2)
    mafftalg=ttk.Combobox(mafftoptframer1,textvariable=MAFFTalgorithm,state='readonly',values=('Auto','Fast (no iterations)','Compromise (2 iterations)','accurate (1000 iterations)','custom'))
    mafftalg.grid(row=1,column=3)
    ttk.Label(mafftoptframer1,text='Distance measure:').grid(row=1,column=4)
    mafftdis=ttk.Combobox(mafftoptframer1,textvariable=MAFFTdistance,state='readonly',values=('6mers (faster)','Local pairwise alignment (very slow, <200 seq)','Global pairwise alignment (very slow, <200 similar-sized seq)','Genaf local pairwise alignment (very slow, <200 seq with large gaps)'),width=60)
    mafftdis.grid(row=1,column=5)
    ttk.Label(mafftoptframer2,text='Iterations:').grid(row=1,column=3)
    mafftitter=ttk.Entry(mafftoptframer2,textvariable=MAFFTiterations,width=20)
    mafftitter.grid(row=1,column=4)
    ttk.Label(mafftoptframer2,text='Additional MAFFT options \n(e.g. --retree 4 --nofft)').grid(row=1,column=1)
    mafftopt=ttk.Entry(mafftoptframer2,textvariable=MAFFToptions,width=20)
    mafftopt.grid(row=1,column=2)
    muscleopt1=ttk.Label(alignoptframe,text='MUSCLE')
    muscleopt1.grid(row=1,column=3)
    muscleopt1.grid_remove()

namealignment=Text(aviewframe,height=10,width=10)
seqalignment=Text(aviewframe,height=10)
namealignment.grid(column=0,row=0)
seqalignment.grid(column=2,row=0)
namescrollv=ttk.Scrollbar(aviewframe,orient=VERTICAL,command=namealignment.yview)
namescrollh=ttk.Scrollbar(aviewframe,orient=HORIZONTAL,command=namealignment.xview)
namescrollh.grid(column=0,row=1,sticky=(N,S,E,W))
def vert_scroll(*args):
    namealignment.yview(*args)
    seqalignment.yview(*args)

alignscrollv=ttk.Scrollbar(aviewframe,orient=VERTICAL,command=vert_scroll)
alignscrollv.grid(column=3,row=0,sticky=(N,S,E,W))
namealignment.configure(xscrollcommand=namescrollh.set,yscrollcommand=alignscrollv.set,wrap='none')
seqscrollh=ttk.Scrollbar(aviewframe,orient=HORIZONTAL,command=seqalignment.xview)
seqscrollh.grid(column=2,row=1,sticky=(N,S,E,W))
seqalignment.configure(xscrollcommand=seqscrollh.set,wrap='none')
slave={namealignment: seqalignment,seqalignment: namealignment}
for i in (namealignment,seqalignment):
    def down(event):
        root.tk.call('tk::TextSetCursor',slave[event.widget],root.tk.call('tk::TextUpDownLine',event.widget,1))
    i.bind('<Down>',down,add=True)
    def up(event):
        root.tk.call('tk::TextSetCursor', slave[event.widget],root.tk.call('tk::TextUpDownLine', event.widget, -1))
    i.bind('<Up>', up, add=True)
    

tbmethod=StringVar()
tbmethod.set('RAxML')
raxmlmodel=StringVar()
raxmlmodel.set('PROTGAMMALG')
raxmodchoose=ttk.Combobox(treebuildframe,textvariable=raxmlmodel,state='readonly',values=('PROTGAMMALG','GTRGAMMA'))
raxmodchoose.grid(row=1,column=2)
tbmchoose=ttk.Combobox(treebuildframe,textvariable=tbmethod,state='readonly',values=('UPGMA','RAxML','Parsimony'))
tbmchoose.grid(row=1,column=1)
tbexe=ttk.Button(treebuildframe,text="Run tree search",command=runtreesearch)
tbexe.grid(row=1,column=3)
#
#
#
#location = ttk.Label(mainframe,textvariable=file_loc)
#location.grid(column=2, row=1, sticky=(W, E))
#
#ttk.Label(mainframe, textvariable=rc).grid(column=2, row=2, sticky=(W, E))
#filebutton=ttk.Button(mainframe,text='openfile', command=filedialog)
#filebutton.grid(column=3, row=1, sticky=(N,W))
#
#ttk.Label(mainframe, text="File location:").grid(column=1, row=1, sticky=W)
#ttk.Label(mainframe, text="Fasta:").grid(column=1, row=2, sticky=E)

for child in mainframe.winfo_children():
    child.grid_configure(padx=5, pady=5)
    child.configure(borderwidth=2,relief='sunken')


#filebutton.focus()
#root.bind('<Return>', inerfunc)

root.mainloop()
