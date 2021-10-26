#-*- coding: utf-8 -*-
"""
Created on Mon. September 24 18:19:38 2020

@author: Song Liu, Tianzhou Ma

"""
######These are packages not use. we may delte these sentences later.
# import os
# import sys
# # 图像处理/展现的相关函数库
# import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
# import matplotlib.colors as colors
# import numpy as np
#for palette as 10 colors.
#import seaborn as sns

from PIL import Image
# Utilities 相关函式库

# 把一些警告的讯息暂时关掉
import warnings
warnings.filterwarnings('ignore')

#process images
from PIL import ImageDraw,ImageFont

#file and download  from web url
import os
import urllib.request

#xml read
from xml.dom.minidom import parse
import xml.dom.minidom

import pandas as pd
from pandas.core.frame import DataFrame

import ast


'''
1) Read from sbgn file to get all genes . 
<glyph class="macromolecule" id="entityVertex_8951627_701">
    <label text="AcM-UBE2M"/>
    <bbox w="56.0" h="38.0" x="571.0" y="225.0"/>
    <glyph class="unit of information" id="entityVertex_8951628_1_mt">
	<label text="Ac"/>
	<bbox w="19.0" h="16.0" x="589.0" y="257.0"/>
    </glyph>
</glyph>
2) Read from file to get the genes that need to be painted
3) Take the intersection of the two lists to decide the genes that can be painted.
4) get the human and mouse color value from the csv to get the corresponding values of the selected genes.
5) According to thhe sbgn elements to get the properties of the gene element to decide the postition, width,height, subunits 
loop the following steps:
6) open the jpeg, use the two colors (half to half) to replace the color of the rectangle area
7) paint the green edges for the recangle
8)add the gene names and the units(copy paste or repaint the unit rectangle) to the rectangle.
end of loop
9) Put a palette image below the image
'''
class ImgText:
  #font = ImageFont.truetype('arial.ttf',8)
  def __init__(self, text,width):
    # 预设宽度 可以修改成你需要的图片宽度
    self.width = width
    # 文本
    self.text = text
    # 段落 , 行数, 行高
    self.duanluo, self.note_height, self.line_height = self.split_text()
  def get_duanluo(self, text):
    txt = Image.new('RGBA', (100, 100), (255, 255, 255, 0))
    draw = ImageDraw.Draw(txt)
    # 所有文字的段落
    duanluo = ""
    # 宽度总和
    sum_width = 0
    # 几行
    line_count = 1
    # 行高
    line_height = 0
    for char in text:
      width, height = draw.textsize(char)
      sum_width += width
      if sum_width > self.width: # 超过预设宽度就修改段落 以及当前行数
        line_count += 1
        sum_width = 0
        duanluo += '\n'
      duanluo += char
      line_height = max(height, line_height)
    if not duanluo.endswith('\n'):
      duanluo += '\n'
    return duanluo, line_height, line_count
  def split_text(self):
    # 按规定宽度分组
    max_line_height, total_lines = 0, 0
    allText = []
    for text in self.text.split('\n'):
      duanluo, line_height, line_count = self.get_duanluo(text)
      max_line_height = max(line_height, max_line_height)
      total_lines += line_count
      allText.append((duanluo, line_count))
    line_height = max_line_height
    total_height = total_lines * line_height
    return allText, total_height, line_height
  # def draw_text(self):
  #   """
  #   这个函数是测试类的绘图以及文字用的。实际不使用。
  #   :return:
  #   """
  #   note_img = Image.open("001.png").convert("RGBA")
  #   draw = ImageDraw.Draw(note_img)
  #   # 左上角开始
  #   x, y = 0, 0
  #   for duanluo, line_count in self.duanluo:
  #     draw.text((x, y), duanluo, fill=(255, 0, 0), font=ImgText.font)
  #     y += self.line_height * line_count
  #   note_img.save("result.png")

#画圆角矩形，imgPath图片路径,color颜色,x横坐标位置,y纵坐标位置,w width,h height,r 圆角半径
def drawRoundRec(drawObject, color, x, y, w, h, r,ifhalf=0):
    '''Rounds'''
    drawObject.ellipse((x, y, x + r, y + r), fill=color)
    drawObject.ellipse((x, y + h - r, x + r, y + h), fill=color)
    if ifhalf==0:
        drawObject.ellipse((x + w - r, y, x + w, y + r), fill=color)
        drawObject.ellipse((x + w - r, y + h - r, x + w, y + h), fill=color)
    else:
        drawObject.rectangle((x + w - r, y, x + w, y + r), fill=color)
        drawObject.rectangle((x + w - r, y + h - r, x + w, y + h), fill=color)

    '''rec.s'''
    drawObject.rectangle((x + r / 2, y, x + w - (r / 2), y + h), fill=color)
    drawObject.rectangle((x, y + r / 2, x + w, y + h - (r / 2)), fill=color)

def paintGenes_and_save(pathwayFileName,Genes,savename,xmax,ymax):
    # # 根目录路径
    # root_dir = os.getcwd()
    # # 训练/验证用的资料目录
    # data_path = os.path.join(root_dir,'data')
    # # 测试用的图像
    # test_image = os.path.join(data_path,'R-HSA-170834.jpeg')

    # 载入图像
    image = Image.open(pathwayFileName)

    # 存储图像并转换格式(jpg->png)
    #image.save(os.path.join(data_path,'new_image.png'))

    #plt.imshow(image)
    #plt.show()

    # new_image = image.resize((400,400))
    #
    # print("原图像的大小:",image.size)
    # print("新图像的大小:",new_image.size)
    #
    # plt.imshow(new_image)
    # plt.show()

    # 定义要裁剪的边界框坐标
    # x1 = 1086.0
    # y1 = 489.0
    # x2 = 1110
    # y2 = 512
    # bbox = (x1,y1,x2,y2)
    #
    # xg1 = 1053
    # yg1 = 432
    # xg2 = 1156
    # yg2 = 499
    left=14
    top=14

    for row in Genes.values:
        print(row[1], row[2])

        xg0=float(row[3])
        yg0=float(row[4])
        gwidth=float(row[1])
        gheight=float(row[2])

        xg1 = xg0+left
        yg1 = yg0+top
        xg2 = xg1+gwidth
        yg2 = yg1+gheight

        #gbbox = (xg1,yg1,xg2,yg2)
        geneTet=row[0]

        hbcolor=row[8]
        mbcolor=row[9]

        # ######the little unit on the box.
        # xu0=589
        # yu0=257
        # uwidth=19
        # uheight=16
        # xu1 =xu0+left
        # yu1 = yu0+top
        # xu2 = xu1+uwidth-4
        # yu2 = yu1+uheight-1
        # ubbox = (xu1,yu1,xu2,yu2)

        # # 进行裁剪
        # cropped_image = image.crop(ubbox)


        img_draw = ImageDraw.Draw(image) # 在blank_image上绘图

        #这里是将旧的画矩形的代码，改为画圆角矩形。
        #这里，通过将所画矩形缩小1号，xg1 + 1, yg1 + 1, gwidth - 1, gheight - 1,露出露出原来对应位置图像的单框线来。
        #这里添加了新的处理，row[7]==4时，unspecified entity时，画的hb和mb对应的矩形要再缩小一些，xg1+6,yg1+6, gwidth-12, gheight-12，以露出原来对应位置图像的双框线来。
        if row[7]==4:
            drawRoundRec(img_draw, mbcolor, xg1+6,yg1+6, gwidth-12, gheight-12, gwidth/12)
            drawRoundRec(img_draw,hbcolor , xg1+6,yg1+6, gwidth/2-12, gheight-12, gwidth/12,1)
            x, y = xg1 + 8, yg1 + 8
            # drawRoundRec(img_draw, mbcolor, xg1 + 1, yg1 + 1, gwidth - 1, gheight - 1, gwidth / 7)
            # drawRoundRec(img_draw, hbcolor, xg1 + 1, yg1 + 1, gwidth / 2 - 1, gheight - 1, gwidth / 7, 1)
            # x, y = xg1 + 3, yg1 + 3
        else:
            drawRoundRec(img_draw, mbcolor, xg1 + 1, yg1 + 1, gwidth - 1, gheight - 1, gwidth / 7)
            drawRoundRec(img_draw, hbcolor, xg1 + 1, yg1 + 1, gwidth / 2 - 1, gheight - 1, gwidth / 7, 1)
            x, y = xg1 + 3, yg1 + 3

        ####旧代码，画矩形的代码
        # img_draw.rectangle((xg1,yg1, xg2, yg2),outline=None,fill=hbcolor)
        # img_draw.rectangle((xg1, yg1, xg1+gwidth/2, yg2),outline=None,fill=mbcolor)
        # img_draw.line([(xg1,yg1),(xg2,yg1),(xg2,yg2),(xg1,yg2),(xg1,yg1)],fill='rgb(104, 136, 112)',width=3)


        #font and write text
        #fnt = ImageFont.truetype('arial.ttf',8)# 修改电脑上的字型,字体可自行下载
        #############################
        #这是原来的写一行字的代码，这里要改为使用上边编写的ImgText类，根据字符长度和矩形长度，自动分段，按行写。
        n = ImgText(geneTet,  gwidth)
        for duanluo, line_count in n.duanluo:
            img_draw.text((x, y), duanluo, fill='rgb(0,2,97)')
            y += n.line_height * line_count
        # img_draw.text((xg1+12,yg1+8),geneTet,font=fnt,fill='rgb(0,2,97)')

    imagePallete = Image.open("pallete.jpeg")
    #medium resolution, not high resolution. need to add a  param.
    print(ymax)
    position = (5, int(ymax+top+5))
    # 进行粘贴
    image.paste(imagePallete,position)

    image.save(savename, quality=95, subsampling=0)

    # plt.imshow(image)
    # plt.show()

'''
递归遍历文件夹，得到所有文件名
'''
def dir_list(path):
    allfile = []
    filelist = os.listdir(path)

    for filename in filelist:
        filepath = os.path.join(path, filename)
        if os.path.isdir(filepath):
            dir_list(filepath, allfile)
        else:
            allfile.append(filepath)

    return allfile

'''
保存远程url图片和数据
'''
def download_and_save(url,savename):
    try:
        urlopen=urllib.request.urlopen(url)
        data = urlopen.read()
        urlopen.close()
        fid=open(savename,"wb+")
        fid.write(data)
        print ("download succeed: "+ url)
        fid.close()
    except IOError:
        print ("download failed: "+ url)


def get_all_imageandxml(datadir,pathwayID):
    sbgnFileURL="https://reactome.org/ContentService/exporter/event/"+pathwayID+".sbgn"
    lowResolutionFileURL="https://reactome.org/ContentService/exporter/diagram/"+pathwayID+".jpeg?diagramProfile=Modern&quality=5"
    highResolutionfile="https://reactome.org/ContentService/exporter/diagram/"+pathwayID+".jpeg?diagramProfile=Modern&quality=7"

    sbgnFileName=pathwayID+".sbgn"
    lowResolutionFileName=pathwayID+".jpeg"
    highResolutionfileName=pathwayID+"(1).jpeg"

    col = 2
    row = 3
    Files = [[0] * col for _ in range(row)]
    Files[0][0]=sbgnFileName
    Files[0][1]=sbgnFileURL
    Files[1][0]=lowResolutionFileName
    Files[1][1]=lowResolutionFileURL
    Files[2][0]=highResolutionfileName
    Files[2][1]=highResolutionfile

    for i in range(len(Files)):
        File_name = Files[i][0]
        File_url = Files[i][1]
        if False == os.path.exists(datadir + '/' + File_name):
            download_and_save(File_url,datadir + '/' + File_name)


def get_allGenes_fromXML(sbgnFileName):
    # use minidom parser to open XML
    DOMTree = xml.dom.minidom.parse(sbgnFileName)
    # get root element
    collection = DOMTree.documentElement
    if collection.hasAttribute("xmlns"):
        print("Root element : %s" % collection.getAttribute("xmlns"))

    # get map
    map = collection.getElementsByTagName("map")

    b = map[0]
    print(b.nodeName, b.getAttribute("language"))

    #Here are the list of all items.
    itemlist = b.getElementsByTagName('glyph')

    col = 8
    row = 1
    genes = [[0] * col for _ in range(row)]
    genesName=[]
    i = 0
    for item in itemlist:
        # "association", "complex", "macromolecule", and "unspecified entity"
        #if (item.getAttribute("class") == "macromolecule") or (item.getAttribute("class") == "association") or (item.getAttribute("class") == "complex") or (item.getAttribute("class") == "unspecified entity") :
        if (item.getAttribute("class") == "macromolecule") or (item.getAttribute("class") == "complex") or (item.getAttribute("class") == "unspecified entity") :
            #get and set label genes:0
            label = item.getElementsByTagName("label")
            value = label[0].getAttribute("text")
            print(value)
            genes[i][0] = value
            genesName.append(value)
            #get and set w h x y genes:1 2 3 4
            bbox = item.getElementsByTagName("bbox")
            w = bbox[0].getAttribute("w")
            h = bbox[0].getAttribute("h")
            x = bbox[0].getAttribute("x")
            y = bbox[0].getAttribute("y")
            genes[i][1] = w
            genes[i][2] = h
            genes[i][3] = x
            genes[i][4] = y

            # get and set graphtype genes:7
            graphtype=0
            if (item.getAttribute("class") == "macromolecule"):
                graphtype = 1
            #elif (item.getAttribute("class") == "association"):
                #graphtype = 2
            elif (item.getAttribute("class") == "complex"):
                graphtype = 3
            elif (item.getAttribute("class") == "unspecified entity"):
                graphtype = 4
            genes[i][7] = graphtype

            genes.append([0, 0, 0, 0, 0, 0, 0,0])
            print(genes[i])
            i = i + 1
    genes = genes[:-1]

    #get the max y of the graph elements for paste pallete image later in the new jpeg.
    ymax=0
    xmax=0

    for item in itemlist:
        bbox = item.getElementsByTagName("bbox")
        w = bbox[0].getAttribute("w")
        h = bbox[0].getAttribute("h")
        x = bbox[0].getAttribute("x")
        y = bbox[0].getAttribute("y")
        if float(x)+float(w) > xmax:
                xmax = float(x)+float(w)
        if float(y)+float(h) > ymax:
                ymax = float(y)+float(h)


    return genes,genesName,xmax,ymax

def CallFromR(signPM,ReactomePath,datadir,pathwayID):
    # function parameters
    #csv_hm_name = "signPM.mat.csv"  # "signPM_hb_mb.csv"
    #pathway_Reactome_csv = ""
    #csv_match_name = "Reactome_all_human_pathways_match.csv"  # "Reactome_all_human_pathways_match.csv"

    imageFile = datadir + pathwayID + ".jpeg"
    NewImageFile = datadir + pathwayID + "New.jpeg"

    # ######download the jpegs and sbgn as the pathwayID.
    get_all_imageandxml(datadir, pathwayID)

    ###get the genes from the sbgn, Here we need to modify as the new rules.
    #get all the genes of "class" == "macromolecule" or  "association" or "complex" or "unspecified entity" from sbgn.
    sbgnFileName = datadir + pathwayID + ".sbgn"
    Genes, genesName, xmax, ymax = get_allGenes_fromXML(sbgnFileName)

    ###get the selected Genes from csv.
    ###pathwayid get genes from csv
    #csvnameNow = csv_match_name
    #dfMatch = pd.read_csv(csvnameNow, sep=',', quotechar='"')  # , index_col=0)
    dfMatch = ReactomePath
    # 解决pymysql.err.InternalError: (1054, "Unknown column 'nan' in 'field list'") 的问题
    # dfMatch2 = dfMatch[dfMatch['Genes'].notnull()]
    dfMatch2 = dfMatch.dropna()
    print("dfMatch2")
    print(dfMatch2)

    SelectedGenes = dfMatch2.loc[dfMatch2['ID'] == pathwayID]
    print("SelectedGenes")
    print(SelectedGenes)

    # If Empty DataFrame
    if SelectedGenes.empty:
        print("No " + pathwayID + " Pathway in the Reactome pathways provided!")

    SelectedGenes = SelectedGenes['Genes'].values.tolist()
    SelectedGenes = SelectedGenes[0].split()
    print("SelectedGenes")
    print(SelectedGenes)

    # Here, the rules need to be modified.不仅要找唯一对应的，要把pathway，list里面的gene name都去框label匹配一遍。有一个genename，用hb,mb；
    # 有两个取平均值。矩形里边还是标原来的label；有几个取几个的平均值。这个怎么匹配，怎么记下来。
    # 先按照完整匹配来做，匹配几个算几个。这个是现在关键的算法。
    #(1)过滤genes，将genes中能匹配到"Reactome_all_human_pathways_match" 对应pathwayid的gene列表的矩形过滤出来。并将对应的匹配基因，添加到gene的[6]中。
    #这里要加特殊处理，一碰到ifn，就加上，不再匹配。只加一次" IFN"。下边算hb和mb，也要特殊处理。
    ifAddIFN=0
    for item in Genes:
        label = item[0]
        item[6] = ""
        for gene in SelectedGenes:
            #" IFN"特殊处理，确保出现的IFN只添加一次。
            if ("IFN" in gene) and ("IFN" in label) and (ifAddIFN == 0):
                ifAddIFN = 1
                item[5] = 1
                item[6] = item[6] + " IFN"
                continue
            elif (gene in label) and not ("IFN" in gene):
                item[5] = 1
                item[6] = item[6] + " " + gene
        ifAddIFN = 0
    df4 = DataFrame(Genes, columns=['label', 'w', 'h', 'x', 'y','ifselected', 'matchgenes', 'graphtype'])
    print("df4")
    print(df4)

    df4 = df4.loc[df4['ifselected'] == 1]
    df4["index"] = range(len(df4))
    df4 = df4.set_index(["index"])
    df4["dat1"] = 0
    df4["dat2"] = 0
    print("df4")
    print(df4)

    ###get the hb and mb values from csv
    #csvnameNow = csv_hm_name
    #df = pd.read_csv(csvnameNow, sep=',', quotechar='"')  # , index_col=0)
    df = signPM
    # 解决pymysql.err.InternalError: (1054, "Unknown column 'nan' in 'field list'") 的问题
    df2 = df.astype(object).where(pd.notnull(df), None)
    print("df2")
    print(df2)
    df2list = df2['Genes'].tolist()

    #计算signPM.mat.csv中IFN系列基因的平均值
    df2hb=0
    df2mb=0
    i = 0
    for item1 in df2.values:
        if ("IFN" in item1[0]):
            i = i + 1
            df2hb =df2hb + item1[1]
            df2mb = df2mb + item1[2]
    if i > 0:
        df2hb = df2hb / i
        df2mb = df2mb / i
    print("df2hb,df2mb: "+str(df2hb)+","+str(df2mb) )

    #(2)过滤genes，将genes中能匹配到"signPM.mat.csv" 对应genes列表的矩形过滤出来。
    for j in range(len(df4)):
        # print(df4.loc[j]['matchgenes'].split())
        # "IFN"特殊处理，先去掉"IFN"，去对应df2匹配子集，同时将"IFN"作为或条件，作为是否选中作为要进行重画的矩形。
        splitstr=df4.loc[j]['matchgenes'].replace(" IFN", " ")
        SelectedGenes = set(splitstr.split())
        if SelectedGenes.issubset(df2list) or (" IFN" in df4.loc[j]['matchgenes']):
            print(SelectedGenes)
            df4.loc[j, 'ifselected'] = 2
    df4 = df4.loc[df4['ifselected'] == 2]
    df4["index"] = range(len(df4))
    df4 = df4.set_index(["index"])
    print("df4")
    print(df4)

    #求genes对应矩形，对应[6]，也就是genes.loc[j]['matchgenes']对应"signPM.mat.csv" 对应gene的hb和mb的平均值。
    #这里，现在要添加一个特殊处理，一碰到ifn，就把表里的signPM.mat.csv所有ifn项都平均。如果还有其他的就再取平均。
    for j in range(len(df4)):
        matchgenes = df4.loc[j]['matchgenes'].split()
        # 去df2中，匹配，匹配一个把hb,mb添加到列表中，最后算平均值
        i = 0
        for item in matchgenes:
            #"IFN"特殊处理，将IFN，对应的df2hb，df2mb加到和中。
            if "IFN" ==item:
                i = i + 1
                df4.loc[j, "dat1"] = df4.loc[j, "dat1"] +df2hb
                df4.loc[j, "dat2"] = df4.loc[j, "dat2"] + df2mb
                continue
            else: #其余的非"IFN"的正常添加对应hb和mb值。
                for item1 in df2.values:
                    if (item == item1[0]):
                        i = i + 1
                        df4.loc[j, "dat1"] = df4.loc[j, "dat1"] + item1[1]
                        df4.loc[j, "dat2"] = df4.loc[j, "dat2"] + item1[2]
        if i > 0:
            df4.loc[j, "dat1"] = df4.loc[j, "dat1"] / i
            df4.loc[j, "dat2"] = df4.loc[j, "dat2"] / i
    print("df4")
    print(df4)

    # #get the corresponding color of hb and mb values from the pallete
    c = ["#1BFC06", "#38EB27", "#53DB45", "#71CA67", "#9CB296", "#AEA7AA", "#BD8688", "#D35556", "#E42D2E",
         "#F70202"]
    # sns.palplot(sns.color_palette(c))
    #
    intervals = [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]

    for j in range(len(df4)):
        print(df4.loc[j, "dat1"], df4.loc[j, "dat2"])
        for i in range(10):
            if intervals[i] <= df4.loc[j, "dat1"] < intervals[i + 1]:
                df4.loc[j, "dat1"] = c[i]
                break
        for i in range(10):
            if intervals[i] <= df4.loc[j, "dat2"] < intervals[i + 1]:
                df4.loc[j, "dat2"] = c[i]
                break
        if df4.loc[j, "dat1"] == 1:
            df4.loc[j, "dat1"] = c[9]
        if df4.loc[j, "dat2"] == 1:
            df4.loc[j, "dat2"] = c[9]
    print("df4")
    print(df4)

    # call the paint function
    paintGenes_and_save(imageFile, df4, NewImageFile, xmax, ymax)


#if __name__ == "__main__":
    #pathwayID = "R-HSA-933541"  # R-HSA-168928 R-HSA-389359    R-HSA-202427   R-HSA-199991
    #datadir = 'images' + '/'
    #CallFromR(datadir,pathwayID)