/*
  Copyright 2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
                 Alex Rand <arand@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeRover.

  VolumeRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifdef USING_TILING

#include <cvc_compat.h>
#include <VolumeRover2/ContourTilerDialog.h>
#include <ContourTiler/tiler.h>
#include <ContourTiler/cl_options.h>

#include <cvcraw_geometry/contours.h>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <QFileDialog>
#include <QMessageBox>

#include "ui_ContourTilerDialog.h"

#include <iostream> 
using namespace std;


//
// 05/01/2011 -- arand -- initial implementation following anisotropic diffusion as a template
// 10/08/2011 -- transfix -- Defaulting to zoomed_volume, and adding extensions from VolMagick

ContourTilerDialog::ContourTilerDialog(QWidget *parent,Qt::WindowFlags flags) 
  : QDialog(parent, flags)
{
  log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeRover2.ContourTilerDialog");

  int idx = -1;

  _ui = new Ui::ContourTilerDialog;
  _ui->setupUi(this);
  
  std::vector<std::string> geoms = 
    cvcapp.data<cvcraw_geometry::contours_t>();
  if (geoms.empty()) {
      QMessageBox::warning(this, tr("Contour Tiler"),
			   tr("No contours loaded."), QMessageBox::Ok);
      _ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);
  }
  else {
    for (const auto& key : geoms)
      _ui->geometryList->addItem(QString::fromStdString(key));

    LOG4CPLUS_TRACE(logger, "found contours: " << geoms.size());

    cvcraw_geometry::contours_t contours =
      cvcapp.data<cvcraw_geometry::contours_t>(geoms[0]);
    LOG4CPLUS_TRACE(logger, "Got something!");
    for (const auto& component : contours.components()) {
      LOG4CPLUS_TRACE(logger, "Component: " << component);
      _ui->componentList->addItem(QString::fromStdString(component));
    }

    _ui->sliceBegin->setText(QString::fromStdString(boost::lexical_cast<string>(contours.z_first())));
    _ui->sliceEnd->setText(QString::fromStdString(boost::lexical_cast<string>(contours.z_last())));

    _ui->additionalArgs->setText("-C 0.01 -e 1e-15 -o raw");

    connect(_ui->geometryList, SIGNAL(currentIndexChanged(int)), SLOT(geometryChangedSlot(int)));
    connect(_ui->outputDirButton,
	    SIGNAL(clicked()),
	    SLOT(outputDirSlot()));
  }
}

ContourTilerDialog::~ContourTilerDialog() {
  delete _ui;
}

void ContourTilerDialog::outputDirSlot() {
  _ui->outputDir->setText(QFileDialog::getExistingDirectory(this,
						    "Output Directory"));
}

void ContourTilerDialog::geometryChangedSlot(int index) {
  log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeRover2.ContourTilerDialog.geometryChangedSlot");

  _ui->componentList->clear();
  if (index == -1) {
    return;
  }

  cvcraw_geometry::contours_t contours =
    cvcapp.data<cvcraw_geometry::contours_t>(_ui->geometryList->itemText(index).toStdString());
  LOG4CPLUS_TRACE(logger, "Got something!");
  for (const auto& component : contours.components()) {
    LOG4CPLUS_TRACE(logger, "Component: " << component);
    _ui->componentList->addItem(QString::fromStdString(component));
  }

  _ui->sliceBegin->setText(QString::fromStdString(boost::lexical_cast<string>(contours.z_first())));
  _ui->sliceEnd->setText(QString::fromStdString(boost::lexical_cast<string>(contours.z_last())));
}

class ContourTilerThread
{
public:
  ContourTilerThread(const vector<Slice>& slices, const Tiler_options& options) : _slices(slices), _options(options) {
    log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeRover2.ContourTilerDialog.ContourTilerThread");
    LOG4CPLUS_TRACE(logger, "Constructor");
  }
  // ContourTilerThread(const cvcraw_geometry::contours_t& contours, const set<string>& components, int first, int last, string args)
  //   : _contours(contours), _components(components), _first(first), _last(last), _args(args) {}

  void operator()()
  {
    log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeRover2.ContourTilerDialog.ContourTilerThread");
    LOG4CPLUS_TRACE(logger, "Execute");

    CVC::ThreadFeedback feedback;

    boost::unordered_map<string, Color> comp2color;
    LOG4CPLUS_TRACE(logger, "Calling tiling");
    tile(_slices.begin(), _slices.end(), comp2color, _options);
    LOG4CPLUS_TRACE(logger, "Tiling complete");
    // if (_output.empty()) {
    //   _output = _volSelected+"_bilateral";
    // }
    // if (_dataset.empty()) {
    //   _dataset = _volSelected+"_bilateral";
    // }


    // // read in the data if necessary
    // if(cvcapp.isData<VolMagick::VolumeFileInfo>(_volSelected))
    //   {
    //     VolMagick::VolumeFileInfo vfi = cvcapp.data<VolMagick::VolumeFileInfo>(_volSelected);

    //     if (_currentIndex == 0) {


    //       if (_output.substr(_output.find_last_of(".")) != _fileType) {
    //         _output = _output + _fileType;   
    //       }

    //       VolMagick::createVolumeFile(cvcapp, _output,
    //                                   vfi.boundingBox(),
    //                                   vfi.dimension(),
    //                                   vfi.voxelTypes(),
    //                                   vfi.numVariables(),
    //                                   vfi.numTimesteps(),
    //                                   vfi.TMin(),
    //                                   vfi.TMax());
    //     }
      
    //     // run anisotropic diffusion
    //     for(unsigned int var=0; var<vfi.numVariables(); var++) {
    //       for(unsigned int time=0; time<vfi.numTimesteps(); time++) {
    //         VolMagick::Volume vol;
	  
    //         readVolumeFile(cvcapp, vol,vfi.filename(),var,time);
    //         vol.bilateralFilter(_radSig,_spatSig,_filRad);
	  
    //         if (_currentIndex == 0) {
    //           writeVolumeFile(cvcapp, vol,_output,var,time);
    //         } else if (_currentIndex == 1) {
    //           // put the dataset in the list
    //           cvcapp.data(_dataset,vol);		  
    //         }
    //       }
    //     }
    //   }
    // else if(cvcapp.isData<VolMagick::Volume>(_volSelected))
    //   {
    //     VolMagick::Volume vol = cvcapp.data<VolMagick::Volume>(_volSelected);
    // 	vol.bilateralFilter(_radSig,_spatSig,_filRad);
      
    //     if (_currentIndex == 0) {
    //       //if _output not set, overwrite the input volume data
    //       if(_output.empty())
    //         cvcapp.data(_volSelected, vol);
    //       else {
    //         _output = _output + _fileType;
    //         cvcapp.data(_output, vol);
    //       }
    //     } else if (_currentIndex == 1) {
    // 	  cvcapp.data(_dataset,vol);
    //     }      
    //   }
  }

private:
  vector<Slice> _slices;
  // cvcraw_geometry::contours_t _contours;
  // set<string> _components;
  // int _first;
  // int _last;
  Tiler_options _options;
};

void ContourTilerDialog::RunContourTiler() {
  // get parameters
  log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeRover2.ContourTilerDialog.RunContourTiler");
  LOG4CPLUS_TRACE(logger, "RUN CONTOUR TILER");

  string dataset = _ui->geometryList->currentText().toStdString();
  cvcraw_geometry::contours_t contours =
    cvcapp.data<cvcraw_geometry::contours_t>(dataset);

  set<string> components;
  QList<QListWidgetItem *> selected = _ui->componentList->selectedItems();
  for (const auto& item : selected) {
    LOG4CPLUS_TRACE(logger, "Component: " << item->text().toStdString());
    components.insert(item->text().toStdString());
  }

  int first = _ui->sliceBegin->displayText().toInt();
  int last = _ui->sliceEnd->displayText().toInt();
  string args = _ui->additionalArgs->displayText().toStdString();

  const std::vector<CONTOURTILER_NAMESPACE::Slice>& all_slices = contours.slices();
  vector<CONTOURTILER_NAMESPACE::Slice> slices;
  for (int i = first; i <= last; ++i) {
    LOG4CPLUS_TRACE(logger, "Adding slice " << i);
    Slice s = all_slices[i-contours.z_first()];
    vector<string> sc;
    s.components(back_inserter(sc));
    if (!components.empty()) {
      for (const auto& c : sc) {
	if (components.find(c) == components.end()) {
	  s.erase(c);
	}
	else {
	  LOG4CPLUS_TRACE(logger, "  Including component: " << c);
	}
      }
    }
    slices.push_back(s);
  }

  boost::char_separator<char> sep(" \t\n\r");
  boost::tokenizer<boost::char_separator<char> > tok(args, sep);
  vector<string> arg_arr;
  arg_arr.push_back("CommandName");
  arg_arr.insert(arg_arr.end(), tok.begin(), tok.end());
  arg_arr.push_back("Filename");
  LOG4CPLUS_TRACE(logger, "Extra arguments:");
  for (const auto& s : arg_arr) {
    LOG4CPLUS_TRACE(logger, "  " << s);
  }
  CONTOURTILER_NAMESPACE::cl_options clo = cl_parse(arg_arr);
  CONTOURTILER_NAMESPACE::Tiler_options options = clo.options;
  // options.output_dir() = "/tmp";
  options.output_dir() = _ui->outputDir->displayText().toStdString();
  options.z_scale() = contours.z_scale();
  LOG4CPLUS_TRACE(logger, "Options:");
  LOG4CPLUS_TRACE(logger, "  out dir: " << options.output_dir());
  LOG4CPLUS_TRACE(logger, "  out raw: " << options.output_raw());

  cvcapp.startThread("contour_tiler_thread_" + dataset,
                     ContourTilerThread(slices, options));

  // double spatSig = _ui->SpatSigEdit->displayText().toDouble();
  // double filRad = _ui->FilRadEdit->displayText().toDouble();

  // // find the volume files to work with
  // std::string volSelected = _ui->VolumeList->currentText().toStdString();
 
  // std::string output = _ui->OutputFilename->displayText().toStdString();
  // std::string dataset = _ui->DataSetName->displayText().toStdString();     

  // cvcapp.startThread("bilateral_filter_thread_" + volSelected,
  //                    ContourTilerThread(radSig, spatSig, filRad,
  //                                               volSelected,
  //                                               output,
  //                                               dataset,
  //                                               _ui->tabWidget->currentIndex(),
                                                // _ui->FileTypeComboBox->currentText().toStdString()));
}


#endif
