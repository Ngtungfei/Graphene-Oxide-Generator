
#import subprocess
import matplotlib.pyplot as plt
from numpy import sin, pi
import random
from copy import deepcopy
#import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import FileDialog
import tkMessageBox

from Tkinter import *
import tkFileDialog

class Example(Frame):

    def __init__(self, parent):
        Frame.__init__(self, parent)   
         
        self.parent = parent
        self.initUI()
        
    def initUI(self):
        
        self.parent.title("GO")
        sw = self.parent.winfo_screenwidth()
        sh = self.parent.winfo_screenheight()
        self.w = sw*0.2
        self.h= 720
        x = int((sw - self.w)/5*4)
        y = (sh - self.h)/2
        self.parent.geometry('%dx%d+%d+%d' % (self.w, self.h, x, y))
        self.pack(fill=BOTH, expand=1)
        
        self.xend=0
        self.yend=0
        self.xst=0
        self.yst=0
        
        
        #self.frame1=Frame(bd=1, relief=SUNKEN)
        #self.frame1.pack(fill=Y,side=LEFT)
        #self.frame2=Frame(bd=1, relief=SUNKEN)
        #self.frame2.pack(fill=Y,side=LEFT)
        #self.frame3=Frame(bd=1, relief=SUNKEN)
        #self.frame3.pack(side=LEFT)
        
        
        self.wb = Label(self, text='1. Generate Graphene:',width=18,anchor=W)
        self.wb.grid(row=0, column=0)
        self.wb = Label(self, text='x (nm)',width=18,anchor=E,fg='red')
        self.wb.grid(row=1, column=0)
        self.wb = Label(self, text='y (nm)',width=18,anchor=E,fg='red')
        self.wb.grid(row=2, column=0)

        self.ww = Scale(self, from_=1, to=10, length=sw*0.1,resolution=0.5, orient=HORIZONTAL, fg='red')
        self.ww.grid(row=1, column=1)
       
        self.www = Scale(self,from_=1, to=10,length=sw*0.1,resolution=0.5, orient=HORIZONTAL, fg='red')
        self.www.grid(row=2, column=1)
        
        #self.quitButton1 = Button(self, text = 'Open xyz.file',width=15,command = self.openwindows)
        #self.quitButton1.grid(row=3, column=0)
        self.quitButton1 = Button(self, text = 'Generate',width=15,command = self.newwindows)
        self.quitButton1.grid(row=3, column=1)
        
        self.wa = Label(self, text='2. Remove Carbons:',width=18,anchor=W)
        self.wa.grid(row=4, column=0)
        self.wb = Label(self, text='Edge C possibility (%)',width=18,anchor=E,fg='red')
        self.wb.grid(row=5, column=0)
        self.wb = Label(self, text='Bulk C possibility (%)',width=18,anchor=E,fg='red')
        self.wb.grid(row=6, column=0)
        self.wb = Label(self, text='Total removing rate (%)',width=18,anchor=E,fg='red')
        self.wb.grid(row=7, column=0)

        self.w1 = Scale(self, from_=0, to=100,length=sw*0.1,orient=HORIZONTAL, fg='red')
        self.w1.grid(row=5, column=1)

        self.w2 = Scale(self, from_=0, to=100,length=sw*0.1,orient=HORIZONTAL, fg='red')
        self.w2.grid(row=6, column=1)

        self.w3 = Scale(self, from_=0, to=100,length=sw*0.1,orient=HORIZONTAL, fg='red')
        self.w3.grid(row=7, column=1)

        self.quitButton2 = Button(self, text = 'OK',width=15,command = self.RemC)
        self.quitButton2.grid(row=8, column=1)
       
        self.wc = Label(self, text='3. -COOH amount (%):',width=18,anchor=W)
        self.wc.grid(row=9, column=0)
        
        self.w4 = Scale(self, from_=0, to=100,length=sw*0.1,orient=HORIZONTAL, fg='blue')
        self.w4.grid(row=9, column=1)
        self.w4.set(100)
        
        
        self.quitButton3 = Button(self, text = 'Add -COOH to Edge',width=15,command = self.addcooh)
        self.quitButton3.grid(row=10, column=1)

        self.wc = Label(self, text='4. -C=O :',width=18,anchor=W)
        self.wc.grid(row=11, column=0)
        
        self.quitButton3 = Button(self, text = 'Add -C=O to Edge',width=15,command = self.addco)
        self.quitButton3.grid(row=12, column=1)

        self.wd = Label(self, text='5. -O- amount (%):',width=18,anchor=W)
        self.wd.grid(row=13, column=0)
        
        self.w6 = Scale(self, from_=0, to=100,length=sw*0.1,orient=HORIZONTAL, fg='blue')
        self.w6.grid(row=13, column=1)
        
        self.quitButton3 = Button(self, text = 'Add -O-',width=15,command = self.addo)
        self.quitButton3.grid(row=14, column=1)

        self.wd = Label(self, text='6. -OH amount (%):',width=18,anchor=W)
        self.wd.grid(row=15, column=0)

        self.w5 = Scale(self, from_=0, to=100,length=sw*0.1,orient=HORIZONTAL, fg='blue')
        self.w5.grid(row=15, column=1)
        
        self.quitButton3 = Button(self, text = 'Add -OH',width=15,command = self.addoh)
        self.quitButton3.grid(row=16, column=1)
        
        self.b1=Button(self, text="Save as",width=15,command= self.savexyz)
        self.b1.grid(row=16, column=0)

        self.CdCbond=1.34#Angstrom
        self.CCbond= 1.54
        self.CdObond=1.23
        self.CObond=1.43
        self.OHbond=0.96
        self.oxide='Graphene'

    def addo(self):
        if self.w6.get()==0:
            return
        self.step=5
        self.delbulkoa=deepcopy(self.delbulkco)
        self.delbulkob=deepcopy(self.delbulkco)
        self.delbulkoha=[]
        self.delbulkohb=[]
        
        self.atomo=deepcopy(self.atomco)#[atom_name x y z atom_type subst_id  subst_name  charge]
        self.atomoh=[]
        
        self.obonding=deepcopy(self.cobonding)#[origin_atom_id target_atom_id bond_type]
        self.ohbonding=[]

        self.addbond=[]
        #self.delbulk=[]
        self.ccbonglist=[]
        for i in self.bonding:
            if i[2]=='2':
                self.ccbonglist.append(i)
                 
        #self.ccbonglist=random.sample(self.bonding,self.ccbongN)
        
        for i in self.ccbonglist:
            if i[0] not in self.delbulkoa and i[1] not in self.delbulkoa:
                self.addbond.append(i)
                self.delbulkoa.append(i[0])
                self.delbulkoa.append(i[1])
       
        n= int(self.w6.get()*len(self.addbond)/100)
        self.olist=random.sample(self.addbond,n)
        
        ab=0
        ca=len(self.atomco)
        
        for i in self.olist:
            self.delbulkob.append(i[0])
            self.delbulkob.append(i[1])
            
            a=float(random.sample([-1.0,1.0],1)[0])
            self.atomo.append(['O',(self.atomc[i[0]][1]+self.atomc[i[1]][1])/2,(self.atomc[i[0]][2]+self.atomc[i[1]][2])/2,a*((self.CObond)**2-(self.CdCbond/2)**2)**0.5,'O.3',1,'UNL1','0.00'])
            self.atomo[i[0]][4]='C.3'
            self.atomo[i[1]][4]='C.3'
            self.obonding.append([i[0], ca+ab,1])
            self.obonding.append([i[1], ca+ab,1]) 
            for ya in self.bonding:
                if ya[0] in i or ya[1] in i:
                    yu =self.bonding.index(ya)
                    #print i, ya,yu
                    self.obonding[yu][2]='1'
                     
            ab+=1
            #upgrade bondlist anglelist
        #print self.bonding
        #print self.obonding
        #print self.atomcooh
        #print self.atomo
         #,self.newbulk
        
        #self.plotall(self.atomo)
        print len(self.olist),'-O- added'
        self.oxide= 'C%dCOOH%dC=O%d-O-%d'%(self.newcarb,self.coohn,len(self.C3),len(self.olist))
        self.plotallbond(self.obonding,self.atomo)
        
        
    def addoh(self):
        if self.w5.get()==0:
            return
        self.step=6
        #self.ohbulk=[]
        self.addbondoh=[]
        self.delbulkoha=deepcopy(self.delbulkob)
        self.delbulkohb=deepcopy(self.delbulkob)
        
        self.atomoh=deepcopy(self.atomo)#[atom_name x y z atom_type subst_id  subst_name  charge]
        
        self.ohbonding=deepcopy(self.obonding)#[origin_atom_id target_atom_id bond_type]
        
        #for i in self.addbond:
            #if i not in self.olist:
               # self.ohbulk.append(i)

        self.ccbonglistoh=[]

        for i in self.bonding:
            if i[2]=='2':
                self.ccbonglistoh.append(i)
        
        
        for i in self.ccbonglistoh:
            if i[0] not in self.delbulkoha and i[1] not in self.delbulkoha:
                self.addbondoh.append(i)
                self.delbulkoha.append(i[0])
                self.delbulkoha.append(i[1])      
               
        n= int(self.w5.get()*len(self.addbondoh)/100)
        self.ohlist=random.sample(self.addbondoh,n)
        
        ca=len(self.atomo)
        ab=0
        
        for i in self.ohlist:
            
            self.delbulkohb.append(i[0])
            self.delbulkohb.append(i[1])
            
            a=float(random.sample([-1.0,1.0],1)[0])
            
            self.atomoh.append(['O',self.atomc[i[0]][1],self.atomc[i[0]][2],a*self.CObond,'O.3',1,'UNL1','0.00'])
            self.ohbonding.append([i[0],ca+ab*4,1])
            
            self.atomoh.append(['H',self.atomc[i[0]][1]+-1*a*self.OHbond,self.atomc[i[0]][2],a*self.CObond,'H',1,'UNL1','0.00'])
            self.ohbonding.append([ca+ab*4,ca+ab*4+1,1])
            
            self.atomoh.append(['O',self.atomc[i[1]][1],self.atomc[i[1]][2],-1*a*self.CObond,'O.3',1,'UNL1','0.00'])
            self.ohbonding.append([i[1],ca+ab*4+2,1])

            self.atomoh.append(['H',self.atomc[i[1]][1]+a*self.OHbond,self.atomc[i[1]][2],-1*a*self.CObond,'H',1,'UNL1','0.00'])
            self.ohbonding.append([ca+ab*4+2,ca+ab*4+3,1])

            self.atomoh[i[0]][4]='C.3'
            self.atomoh[i[1]][4]='C.3'
            for ya in self.bonding:
                if ya[0] in i or ya[1] in i:
                    yu =int(self.bonding.index(ya))
                    #print i, ya,yu
                    self.ohbonding[yu][2]='1'
    
            ab+=1  
            #upgrade bondlist anglelist
        x=(self.coohn+len(self.C3)+len(self.olist)*2+len(self.ohlist)*2)/float(self.newcarb)*100
        print len(self.ohlist)*2,'-OH added' #,self.newbulk
        #self.oxide= ' %.1f%% oxidation'%x
        self.oxide= 'C%dCOOH%dC=O%d-O-%dOH%d-%.1f%% oxidation'%(self.newcarb,self.coohn,len(self.C3),len(self.olist),len(self.ohlist)*2,x)
        print 'Total %.1f '%x,'% carbons oxidized!'
        #self.plotall(self.atomoh)
        self.plotallbond(self.ohbonding,self.atomoh)
        
        
    def addcooh(self):
        if self.w4.get()==0:
            return
        self.step=3

        self.newedge=deepcopy(self.edge)
        
        for i in self.C3:
            self.newedge.remove(i)
        
        self.atomcooh=deepcopy(self.atomc)#[atom_name x y z atom_type subst_id  subst_name  charge]
        self.atomo=[]
        self.atomoh=[]
        self.atomco=[]
        
        self.coohbonding=deepcopy(self.bonding)#[origin_atom_id target_atom_id bond_type]
        self.obonding=[]
        self.ohbonding=[]
        self.cobonding=[]
        
        
        self.delbulk=[]
        self.delbulkoa=[]
        self.delbulkob=[]
        self.delbulkoha=[]
        self.delbulkohb=[]
        self.delbulkco=[]
        
        m=0
        
        n=int(self.w4.get()*len(self.newedge)/100)
        self.coohn=n

        self.coohlist=random.sample(self.newedge,n)#self.coohlist carbon no linked with COOH
        
        #print self.atomc[4],edge,self.coohlist,self.bonding, self.atomc[1],self.atomc[2],self.atomcooh
        ab=0
        ca=len(self.atomc)
        
        for i in self.coohlist:
           x=0
           y=0

           for b in self.bonding:
                if b[0] == i:
                    x+=self.atomc[b[1]][1]/2
                    y+=self.atomc[b[1]][2]/2
                if b[1] == i:
                    x+=self.atomc[b[0]][1]/2
                    y+=self.atomc[b[0]][2]/2


           co=(self.atomc[i][1]-x)/((self.atomc[i][2]-y)**2+(self.atomc[i][1]-x)**2)**0.5
           si=(self.atomc[i][2]-y)/((self.atomc[i][2]-y)**2+(self.atomc[i][1]-x)**2)**0.5
           self.delbulk.append(i)
           

           
           self.atomcooh.append(['C',self.atomc[i][1]+self.CCbond*co,self.atomc[i][2]+self.CCbond*si,0,'C.2',1,'UNL1','0.00'])
           self.coohbonding.append([i, ca+4*ab,1])
           self.atomcooh.append(['O',self.atomc[i][1]+(self.CCbond+self.CdObond/2)*co,self.atomc[i][2]+(self.CCbond+self.CdObond/2)*si,self.CdObond*sin(pi/3),'O.2',1,'UNL1','0.00'])
           self.coohbonding.append([ca+4*ab, ca+4*ab+1,2])
           self.atomcooh.append(['O',self.atomc[i][1]+(self.CCbond+self.CObond/2)*co,self.atomc[i][2]+(self.CCbond+self.CObond/2)*si,-self.CObond*sin(pi/3),'O.3',1,'UNL1','0.00'])
           self.coohbonding.append([ca+4*ab, ca+4*ab+2,1])
           self.atomcooh.append(['H',self.atomc[i][1]+(self.CCbond+self.CObond/2+self.OHbond)*co,self.atomc[i][2]+(self.CCbond+self.CObond/2+self.OHbond)*si,-self.CObond*sin(pi/3),'H',1,'UNL1','0.00'])
           self.coohbonding.append([ca+4*ab+2, ca+4*ab+3,1])
           ab+=1

        self.newcarb+=n  
        #print self.coohbonding
        #self.plotall(self.atomcooh)
        print '%d -COOH added, still %d edge carbons left'%(n,len(self.edge)-n)
        self.oxide= 'C%dCOOH%d'%(self.newcarb,self.coohn)
        self.plotallbond(self.coohbonding,self.atomcooh)
        
    def addco(self):
        
        self.step=4
        self.atomco=deepcopy(self.atomcooh)
        self.cobonding=deepcopy(self.coohbonding)

        
        ab=0
        ca=len(self.atomcooh)
        
        for i in self.C3:
           x=0
           y=0

           for b in self.bonding:
                if b[0] == i:
                    x+=self.atomc[b[1]][1]/2
                    y+=self.atomc[b[1]][2]/2
                if b[1] == i:
                    x+=self.atomc[b[0]][1]/2
                    y+=self.atomc[b[0]][2]/2


           co=(self.atomc[i][1]-x)/((self.atomc[i][2]-y)**2+(self.atomc[i][1]-x)**2)**0.5
           si=(self.atomc[i][2]-y)/((self.atomc[i][2]-y)**2+(self.atomc[i][1]-x)**2)**0.5
           self.delbulkco.append(i)
           
           self.atomco.append(['O',self.atomc[i][1]+self.CdObond*co,self.atomc[i][2]+self.CdObond*si,0,'O.2',1,'UNL1','0.00'])
           self.cobonding.append([i, ca+ab,2])
           self.atomco[i][4]='C.2'
           ab+=1
           for ya in self.bonding:
                if ya[0] == i or ya[1] == i:
                    yu =int(self.bonding.index(ya))
                    #print i, ya,yu
                    self.cobonding[yu][2]='1'
                    
        print '%d -C=O added'%(len(self.C3))
        self.oxide= 'C%dCOOH%dC=O%d'%(self.newcarb,self.coohn,len(self.C3))
        self.plotallbond(self.cobonding,self.atomco)
            
             
        
    def savexyz(self):
        try:
            if self.step ==1:
                ftypes=[('XYZ file','.xyz')]
                f = tkFileDialog.asksaveasfilename(filetypes=ftypes,initialfile='Graphene.xyz')
                file = open(f, "w")
                file.write('%d'%len(self.carbonx))
                file.write('\n    By WTF')
                for i in range (0,len(self.carbonx)):
                     file.write('\nC     '+str(self.carbonx[i])+'      '+str(self.carbony[i])+'      '+str(self.carbonz[i]))
                file.close()

                print 'One xyz file saved!'

            else:
                         
               if self.step ==2:
                   atoms=list(self.atomc)
                   bonds=list(self.bonding)
                   filename='Graphene(%s).mol2'%self.oxide
               if self.step ==3:
                   atoms=list(self.atomcooh)
                   bonds=list(self.coohbonding)
                   filename='GO(%s).mol2'%self.oxide
               if self.step ==4:
                   atoms=list(self.atomco)
                   bonds=list(self.cobonding)
                   filename='GO(%s).mol2'%self.oxide
               if self.step ==5:
                   atoms=list(self.atomo)
                   bonds=list(self.obonding)
                   filename='GO(%s).mol2'%self.oxide
               if self.step ==6:
                   atoms=list(self.atomoh)
                   bonds=list(self.ohbonding)
                   filename='GO(%s).mol2'%self.oxide

               ftypes=[('Tripo mol2 file','.mol2')]
               f = tkFileDialog.asksaveasfilename(filetypes=ftypes,initialfile=filename)
               file = open(f, "w")
        
               file.write('@<TRIPOS>MOLECULE')
               file.write('\nGraphene Oxide')
               file.write('\n %d %d 0 0 0'%(len(atoms),len(bonds)))
               file.write('\nSMALL')
               file.write('\nGASTEIGER')
               file.write('\n')
               file.write('\n@<TRIPOS>ATOM')

               a=1
               for i in atoms:
                   file.write('\n      %d %s      %.4f   %.4f   %.4f %s   %d  %s      %s'%(a,i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]))
                   a+=1
            
               a=1
               file.write('\n@<TRIPOS>BOND')
               for i in bonds:
                   file.write('\n      %d   %d   %d   %s'%(a,i[0]+1,i[1]+1,i[2]))
                   a+=1
        
               file.close()
               print 'One mol2 file saved! Further geometry optimization is necessary!!!'
        except:
                pass

        
    def RemC(self):
      
        if self.w1.get()+self.w2.get()==0:
            return
        if self.w3.get()==0:
            return

        self.step=2
        self.newcarbonbond=deepcopy(self.carbonbond)
        
        self.newbonds=deepcopy(self.bonds)# carbon bonding no left after removing
        
        self.total= []# carbon no left after removing
        
        for i in range (0,self.carbonN):
            self.total.append(i)
        
        m=0
        
        while m*100/self.carbonN <self.w3.get():

            a=random.choice(self.total)
            if self.newcarbonbond[a] == 2 and random.randrange(0,100) <  self.w1.get():
                m=self.remove(a,m)
            if self.newcarbonbond[a] == 3 and random.randrange(0,100) <  self.w2.get():
                m=self.remove(a,m)
                    
            while 1 in self.newcarbonbond:
                a=self.newcarbonbond.index(1)
                m=self.remove(a,m)

        self.bonds.sort()
        for i in self.bonds:
            if i[0] not in self.total or i[1] not in self.total:
                self.newbonds.remove(i)
            else:
                if i[0]>i[1]:
                    self.newbonds.remove(i)        
        
        self.atomc=[]#final [atom_name x y z atom_type subst_id  subst_name  charge]
        self.atomcooh=[]
        self.atomo=[]
        self.atomoh=[]
        self.atomco=[]

        self.atomce=[]#self.atomc[4]
        
        self.bonding=[] # [origin_atom_id target_atom_id bond_type]
        self.coohbonding=[]
        self.obonding=[]
        self.ohbonding=[]
        self.cobonding=[]
        
        
        self.dict={}
        i=0
        for n in self.total:
            self.atomc.append(['C',self.carbonx[n],self.carbony[n],self.carbonz[n],'C.ar',1,'UNL1','0.00'])
            self.dict[n]=i
            i+=1
        for n in self.newcarbonbond:
            if n> 1:
                   self.atomce.append(n)
            
        for n in self.newbonds:
            self.bonding.append([self.dict[n[0]],self.dict[n[1]],'ar'])
                        
       
        
        self.ccbongN=len(self.bonding)
        self.edge=[]
        self.bulk=[]
        
        for i in range(0, len(self.atomce)):
            if self.atomce[i]==2:
                self.edge.append(i)
                self.bulk.append(i)
            if self.atomce[i]==3:
                self.bulk.append(i)
                
        #print  self.atomc,self.bonding
        #self.plotall(self.atomc)
        print '%d carbons preserved for adding functional groups in next steps'%len(self.total)
        self.oxide='C'+str(len(self.total))
        #print self.bonding
        self.changeartoc2()
        self.C3=[]

        for i in self.edge:
            a=0
            for hi in self.bonding:
                if i in hi and hi[2]=='1':
                    a+=1
                if i in hi and hi[2]=='2':
                    a-=1
            if a==2:
                self.C3.append(i)
                
        #print self.C3,self.edge
              
        self.newcarb=len(self.total)
        
        self.plotallbond(self.bonding,self.atomc)

    def changeartoc2(self):
        listy=[]
        for i in self.atomc:
            listy.append(i[2])
        yy=list(set(listy))
        yy.sort()
        ygroup=[]

        for i in yy:
            ygroup.append([])
            
        for i in range (0,len(self.atomc)):
            ygroup[yy.index(self.atomc[i][2])].append(i)
                
        for i in ygroup:
            if ygroup.index(i)%3==0:
                for jh in range(0,len(i)-1):
                    if [i[jh],i[jh+1],'ar'] in self.bonding:
                        if self.atomc[i[jh]][4]!='C.2'and self.atomc[i[jh+1]][4]!='C.2':
                            self.bonding[self.bonding.index([i[jh],i[jh+1],'ar'])][2]='2'
                            self.atomc[i[jh]][4]='C.2'
                            self.atomc[i[jh+1]][4]='C.2'
            
                        

            if ygroup.index(i)%3==1:
                try:
                    for jh in i:
                       for jk in ygroup[ygroup.index(i)+1]:
                           if [jh,jk,'ar'] in self.bonding:
                               if self.atomc[jk][4]!='C.2':
                                  self.bonding[self.bonding.index([jh,jk,'ar'])][2]='2'
                                  self.atomc[jh][4]='C.2'
                                  self.atomc[jk][4]='C.2'
                except:
                    pass
                try:
                    for jh in i:
                       for jk in ygroup[ygroup.index(i)-1]:
                           if [jk,jh,'ar'] in self.bonding:
                               if self.atomc[jk][4]!='C.2'and self.atomc[jh][4]!='C.2':
                                  self.bonding[self.bonding.index([jk,jh,'ar'])][2]='2'
                                  self.atomc[jh][4]='C.2'
                                  self.atomc[jk][4]='C.2'
                except:
                    pass

                try:
                  for jh in i:

                       if self.atomc[jh][4]=='C.ar'and self.atomc[i[i.index(jh)+1]][4]=='C.ar':
                             if [jh,i[i.index(jh)+1],'ar'] in self.bonding:
                                          self.bonding[self.bonding.index([jh,i[i.index(jh)+1],'ar'])][2]='2'
                                          self.atomc[jh][4]='C.2'
                                          self.atomc[i[i.index(jh)+1]][4]='C.2'
                except:
                    pass

                
            if ygroup.index(i)%3==2:  
                if self.atomc[i[len(i)-1]][4]=='C.ar':
                    try:
                        for jk in ygroup[ygroup.index(i)+1]:
                                if [i[len(i)-1],jk,'ar'] in self.bonding:
                                    if self.atomc[jk][4]!='C.2':
                                          self.bonding[self.bonding.index([i[len(i)-1],jk,'ar'])][2]='2'
                                          self.atomc[i[len(i)-1]][4]='C.2'
                                          self.atomc[jk][4]='C.2'
                    except:
                        pass
                    try:
                        for jk in ygroup[ygroup.index(i)-1]:
                                if [jk,i[len(i)-1],'ar'] in self.bonding:
                                    if self.atomc[jk][4]!='C.2':
                                          self.bonding[self.bonding.index([jk,i[len(i)-1],'ar'])][2]='2'
                                          self.atomc[i[len(i)-1]][4]='C.2'
                                          self.atomc[jk][4]='C.2'
                    except:
                        pass

                if self.atomc[i[0]][4]=='C.ar':
                    try:
                        for jk in ygroup[ygroup.index(i)+1]:
                                 if [i[0],jk,'ar'] in self.bonding:
                                     if self.atomc[jk][4]!='C.2':
                                          self.bonding[self.bonding.index([i[0],jk,'ar'])][2]='2'
                                          self.atomc[i[0]][4]='C.2'
                                          self.atomc[jk][4]='C.2'
                    except:
                            pass
                    try:
                        for jk in ygroup[ygroup.index(i)-1]:
                                 if [jk,i[0],'ar'] in self.bonding:
                                     if self.atomc[jk][4]!='C.2':
                                           self.bonding[self.bonding.index([jk,i[0],'ar'])][2]='2'
                                           self.atomc[i[0]][4]='C.2'
                                           self.atomc[jk][4]='C.2'
                    except:
                            pass

                if self.atomc[i[len(i)-1]][4]=='C.ar'and self.atomc[i[len(i)-2]][4]=='C.ar':
                    if [i[len(i)-2],i[len(i)-1],'ar'] in self.bonding:
                                          self.bonding[self.bonding.index([i[len(i)-2],i[len(i)-1],'ar'])][2]='2'
                                          self.atomc[i[len(i)-1]][4]='C.2'
                                          self.atomc[i[len(i)-2]][4]='C.2'
                if self.atomc[i[0]][4]=='C.ar'and self.atomc[i[1]][4]=='C.ar':
                    if [i[0],i[1],'ar'] in self.bonding:
                                          self.bonding[self.bonding.index([i[0],i[1],'ar'])][2]='2'
                                          self.atomc[i[0]][4]='C.2'
                                          self.atomc[i[2]][4]='C.2'
                                 
            if ygroup.index(i)==len(i)-1:
                for jh in range(0,len(i)-1):
                    if [i[jh],i[jh+1],'ar'] in self.bonding:
                        if self.atomc[i[jh]][4]!='C.2'and self.atomc[i[jh+1]][4]!='C.2':
                            self.bonding[self.bonding.index([i[jh],i[jh+1],'ar'])][2]='2'
                            self.atomc[i[jh]][4]='C.2'
                            self.atomc[i[jh+1]][4]='C.2'
                            
        for i in self.bonding:
             if i[2]=='ar':
                   i[2]='1'

        for i in self.atomc:
             if i[4]=='C.ar':
                   i[4]='C.2'

        #print len(ygroup),self.bonding
            
        
        
        
        
    def plotall(self,li):
        col={'C':'r','O':'b','H':'g'}
        siz={'C':'600','O':'650','H':'300'}
        try:
            self.fig.clf()
        except:
            self.fig = plt.figure(1)

       
        for i in li:
                plt.scatter(i[1], i[2],color=col[i[0]],alpha=.4,s=int(siz[i[0]])/self.ww.get()/self.www.get())

        plt.xlim([self.xst,self.xend])
        plt.ylim([self.yst,self.yend])
        plt.axes().set_aspect('equal')
        plt.xlabel('Angstrom')
        plt.ylabel('Angstrom')
        plt.title('Carbon Number: '+str(len(self.total)))
        plt.show()

    def plotallbond(self,bli,ali):
        col={'C':'r','O':'b','H':'g'}
        siz={'C':'600','O':'650','H':'300'}
        try:
            self.fig.clf()
        except:
            self.fig = plt.figure(1)

       
        for i in bli:
                plt.scatter(ali[i[0]][1], ali[i[0]][2],label=ali[i[0]][0],color=col[ali[i[0]][0]],alpha=.4,s=int(siz[ali[i[0]][0]])/self.ww.get()/self.www.get())
                plt.scatter(ali[i[1]][1], ali[i[1]][2],label=ali[i[1]][0],color=col[ali[i[1]][0]],alpha=.4,s=int(siz[ali[i[1]][0]])/self.ww.get()/self.www.get())
                plt.plot([ali[i[0]][1],ali[i[1]][1]], [ali[i[0]][2],ali[i[1]][2]],color='k')


        plt.xlim([self.xst,self.xend])
        plt.ylim([self.yst,self.yend])

        patchC = Line2D(range(1), range(1), color="white", marker='o', markersize=10,markeredgecolor="red", markerfacecolor="red",alpha=0.8)
        patchO = Line2D(range(1), range(1), color="white", marker='o', markersize=11,markeredgecolor="blue", markerfacecolor="blue",alpha=0.8)
        patchH = Line2D(range(1), range(1), color="white", marker='o', markersize=6, markeredgecolor="green",markerfacecolor="green",alpha=0.4) #mpatches.Patch(color='green',alpha=.4,label='H')
        plt.legend([patchC,patchO,patchH],['C','O','H'],loc='lower left',bbox_to_anchor=(1, 0),numpoints = 1)
        
        plt.axes().set_aspect('equal')
        plt.xlabel('Angstrom')
        plt.ylabel('Angstrom')
        plt.title('Carbon Number: '+str(self.newcarb))
        plt.show()
  
    def remove(self,a,m):
        self.newcarbonbond[a]=0
        m+=1
        self.total.remove(a)
                                 
        if a//(2*self.x)%2==1:
                        f=1
        else:
                        f=-1

        if a%(2*self.x)%2==0 and a-f >=0 and a+f < self.carbonN:
                        if f==1 and a%(2*self.x)==0:
                           pass
                        else:
                           self.newcarbonbond[a-f]-=1
                           if self.newcarbonbond[a-f] ==0:
                                 m+=1
                                 self.total.remove(a-f)
                                          
                           
        if a%(2*self.x)%2==1 and a+f >=0 and a+f < self.carbonN:
                       if f==1 and a%(2*self.x)==2*self.x-1:
                           pass
                       else:
                           self.newcarbonbond[a+f]-=1
                           if self.newcarbonbond[a+f] ==0:
                                 m+=1
                                 self.total.remove(a+f)
                           
            
        if a+2*self.x < self.carbonN:
                        self.newcarbonbond[a+2*self.x]-=1
                        if self.newcarbonbond[a+2*self.x] ==0:
                                m+=1
                                self.total.remove(a+2*self.x)
                                       
                        
        if a-2*self.x>=0:
                        self.newcarbonbond[a-2*self.x]-=1
                        if self.newcarbonbond[a-2*self.x] ==0:
                               m+=1
                               self.total.remove(a-2*self.x)   
        
        return m
    
    def newwindows(self):
        self.step=1
        
        self.carbonx=[]
        self.carbony=[]
        self.carbonz=[]
        self.carbonbond=[]
        self.bonds=[]
        self.carbonN=0
        self.x=int((self.ww.get()*10)//(3*self.CdCbond))
        self.y=int((self.www.get()*10)//(sin(pi/3)*self.CdCbond))
        
        self.atomc=[]#final [atom_name x y z atom_type subst_id  subst_name  charge]
        self.atomcooh=[]
        self.atomo=[]
        self.atomoh=[]

        self.atomce=[]#self.atomc[4]
        
        self.bonding=[] # [origin_atom_id target_atom_id bond_type]
        self.coohbonding=[]
        self.obonding=[]
        self.ohbonding=[]

        if self.y%2==0:
            self.y+=1

        for i in range(0,self.y):
            if i%2==0:
                for m in range(1,self.x+1):
                      self.carbonx.append(m*3*self.CdCbond-self.CdCbond)
                      self.carbony.append(i*sin(pi/3)*self.CdCbond)
                      self.carbonz.append(0)
                      self.carbonN+=1
                      if i==0 or i== self.y-1:
                          self.carbonbond.append(2)
                      else:
                          self.carbonbond.append(3)
                          
                      self.carbonx.append(m*3*self.CdCbond)
                      self.carbony.append(i*sin(pi/3)*self.CdCbond)
                      self.carbonz.append(0)
                      self.carbonN+=1
                      if i==0 or i== self.y-1:
                          self.carbonbond.append(2)
                      else:
                          self.carbonbond.append(3)
        
            else:
                for m in range(1,self.x+1):
                    self.carbonx.append(m*3*self.CdCbond-1.5*self.CdCbond)
                    self.carbony.append(i*sin(pi/3)*self.CdCbond)
                    self.carbonz.append(0)
                    self.carbonN+=1
                    if m==1:
                          self.carbonbond.append(2)
                    else:
                          self.carbonbond.append(3)
                    
                    self.carbonx.append(m*3*self.CdCbond+0.5*self.CdCbond)
                    self.carbony.append(i*sin(pi/3)*self.CdCbond)
                    self.carbonz.append(0)
                    self.carbonN+=1
                    if m==self.x:
                          self.carbonbond.append(2)
                    else:
                          self.carbonbond.append(3)
                      
        self.xend=max(self.carbonx)+3*self.CdCbond
        self.yend=max(self.carbony)+3*self.CdCbond
        self.xst=min(self.carbonx)-3*self.CdCbond
        self.yst=min(self.carbony)-3*self.CdCbond
        
        #print str(self.carbonN), 'Carbon generated'

        for i in range(0,self.carbonN):
             a=i
            
             if a//(2*self.x)%2==1:
                        f=1
             else:
                        f=-1
                           
             if a%(2*self.x)%2==0 and a-f >=0 and a+f < self.carbonN:
                       if f==1 and a%(2*self.x)==0:
                           pass
                       else:
                           #print a-f
                           
                           self.bonds.append([a,a-f])
                                       
             if a%(2*self.x)%2==1 and a+f >=0 and a+f < self.carbonN:
                       if f==1 and a%(2*self.x)==2*self.x-1:
                           pass
                       else:
                           #print a+f
                           self.bonds.append([a,a+f])
                           
             if a+2*self.x < self.carbonN:
                        
                        #print a+2*self.x
                        self.bonds.append([a,a+2*self.x])
                        
             if a-2*self.x>=0:
                        
                        #print a-2*self.x
                        self.bonds.append([a,a-2*self.x])
                        
       
        try:
            self.fig.clf()
        except:
            self.fig = plt.figure(1)
        plt.scatter(self.carbonx, self.carbony,label='C',color='r',alpha=.8,s=600/self.ww.get()/self.www.get())
        plt.xlim([self.xst,self.xend])
        plt.ylim([self.yst,self.yend])
        plt.axes().set_aspect('equal')
        plt.xlabel('Angstrom')
        plt.ylabel('Angstrom')
        plt.title('Carbon Number: '+str(self.carbonN))
        patchC = Line2D(range(1), range(1), color="white", marker='o', markersize=10,markeredgecolor="red", markerfacecolor="red",alpha=0.8)
        patchO = Line2D(range(1), range(1), color="white", marker='o', markersize=11,markeredgecolor="blue", markerfacecolor="blue",alpha=0.8)
        patchH = Line2D(range(1), range(1), color="white", marker='o', markersize=6, markeredgecolor="green",markerfacecolor="green",alpha=0.4) #mpatches.Patch(color='green',alpha=.4,label='H')
        plt.legend([patchC,patchO,patchH],['C','O','H'],loc='lower left',bbox_to_anchor=(1, 0),numpoints = 1)
        print '**!! %d carbons generated for cutting to draw GO'%self.carbonN
        plt.show()

        



def main():
  
    root = Tk()
    ex = Example(root)
    root.mainloop()
    


if __name__ == '__main__':
    main()  



    
    


        
 

