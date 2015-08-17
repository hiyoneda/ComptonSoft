#!/usr/bin/env ruby

require 'ComptonSoftLib'

Degree = Math::PI / 180.0

class MyApp < ANL::ANLApp
  attr_accessor :inputs, :output

  def setup()
    add_namespace ComptonSoft

    chain :CSHitCollection
    chain :DetectorGroupManager
    with_parameters(filename: "database/detector_group.txt")
    chain :ReadComptonEventTree
    with_parameters(file_list: inputs)
    chain :HistogramAzimuthAngle
    with_parameters(number_of_bins: 72)
    chain :SaveData
    with_parameters(output: output)
  end
end

app = MyApp.new
app.inputs = ["compton_cut.root"]
app.output = "compton_phi.root"
app.run(:all, 10000)
