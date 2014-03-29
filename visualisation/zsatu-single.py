#!/usr/bin/python

import vtk
import vtk.util.numpy_support
import numpy
import os
import sys

SRC_DIR = "./x-data/"

def captureImage(window, timestamp):
        image = vtk.vtkWindowToImageFilter()
        image.SetInput(window)
        image.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName("./saturation-{}.png".format(timestamp))
        writer.SetInput(image.GetOutput())

        #window.Render()
        writer.Write()


def addPlot(chart, reader, name):
        data1 = reader.GetOutput()
        
        coords = vtk.vtkFloatArray()
        coords.SetName("Coords")

        for i in range(data1.GetNumberOfPoints()):
                x,y,z = data1.GetPoint(i)
		coords.InsertNextValue(y)

	table = vtk.vtkTable()
	table.AddColumn(coords)
	table.AddColumn(data1.GetPointData().GetArray("S_g"))
        table.AddColumn(data1.GetPointData().GetArray("S_n"))
        table.AddColumn(data1.GetPointData().GetArray("S_w"))

        line1 = chart.AddPlot(vtk.vtkChart.LINE)
        line1.SetInput(table, 0, 1)
	line1.SetMarkerStyle(vtk.vtkPlotPoints.CROSS)

        line2 = chart.AddPlot(vtk.vtkChart.LINE)
        line2.SetInput(table, 0, 2)
	line2.SetMarkerStyle(vtk.vtkPlotPoints.PLUS)

        line3 = chart.AddPlot(vtk.vtkChart.LINE)
        line3.SetInput(table, 0, 3)
        line3.SetMarkerStyle(vtk.vtkPlotPoints.CIRCLE)

        #line1.GetPen().SetLineType(vtk.vtkPen.NO_PEN)
        #line1.SetMarkerStyle(vtk.vtkPlotPoints.PLUS)
        #line1.SetColor(150, 100, 0, 255)

if __name__ == "__main__":
	dirname, filename = os.path.split(sys.argv[1])
        _, _, k = filename.partition("-")
        l,_,_ = k.partition(".")

        view = vtk.vtkContextView()
        view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
        view.GetRenderWindow().SetSize(800, 600)
        chart = vtk.vtkChartXY()
        view.GetScene().AddItem(chart)

	chart.SetTitle("time = {}".format(int(l)))
	chart.GetAxis(0).SetTitle("Saturation")
	chart.GetAxis(1).SetTitle("Coordinate")

	reader = vtk.vtkXMLStructuredGridReader()
	reader.SetFileName(sys.argv[1])
	reader.Update()

	addPlot(chart, reader, "time = {}".format(int(l)))

        chart.SetShowLegend(True)
        view.GetRenderWindow().SetMultiSamples(0)

        view.GetRenderWindow().Render()

        captureImage(view.GetRenderWindow(), l)

