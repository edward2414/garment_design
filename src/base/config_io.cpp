//
//  Filename         : config_io.cpp
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Configuration management class with I/O capabilities.
//  Date of creation : 26/02/2003
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <iostream>
#include <qfileinfo.h>
#include <qdir.h>
#include "config_io.h"

ConfigIO::ConfigIO(QString filename, const QString& doc_type, bool automatic, const QString& sep)
  : _default_file(filename), _automatic(automatic) {
  _doc_type = doc_type;
  _path_sep = sep;
  if (_automatic)
    loadFile();
}

ConfigIO::~ConfigIO() {
  if (_automatic)
    saveFile();
}

QString	ConfigIO::getDefaultFile() const {
  return _default_file;
}

void	ConfigIO::setDefaultFile(const QString& filename) {
  _default_file = filename;
}

bool	ConfigIO::getAuto() const {
  return _automatic;
}

void	ConfigIO::setAuto(bool automatic) {
  _automatic = automatic;
}

QString	ConfigIO::getPathSep() const {
  return _path_sep;
}

void	ConfigIO::setPathSep(const QString& sep) {
  _path_sep = sep;
}

int	ConfigIO::loadFile(const QString& _filename) {

  const QString filename = _filename.isEmpty() ? _default_file : _filename;

  // check wether filename is a valid file and is readable
  QFileInfo fileinfo(filename);
  if (!fileinfo.isFile() || !fileinfo.isReadable()) {
    std::cerr << "Warning: unable to load configuration file \""
	 << fileinfo.fileName() << "\"" << std::endl;
    return 1;
  }

  // read the DOM tree from file
  QFile file(filename);
  file.open(IO_ReadOnly);
  _tree.setContent(&file);
  file.close();

  return 0;
}

int	ConfigIO::saveFile(const QString& _filename) const {

  QString str_tree = _tree.toString();
  if (str_tree.isEmpty())
    return 1;

  const QString filename = _filename.isEmpty() ? _default_file : _filename;

  // if the directory in which we want to generate a file
  // does not exist yet, try to create it
  QFileInfo fileinfo(filename);
  if (!fileinfo.exists()) {
    QDir dir;
    dir.mkdir(fileinfo.dirPath());
  }

  // check wether filename is a valid file and is writable
  QFile file(filename);
  file.open(IO_WriteOnly);
  if (!file.isOpen()) {
    std::cerr << "Warning: unable to save configuration file \""
	 << fileinfo.fileName() << "\"" << std::endl;
    return 1;
  }

  // write the DOM tree to file
  QTextStream out(&file);
  out << str_tree;
  file.close();

  return 0;
}
