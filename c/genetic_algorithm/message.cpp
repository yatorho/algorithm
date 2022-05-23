//
// Created by 32270 on 2022/3/8.
// Author yatorho

#include "message.h"


bool MessageBox::addText(char *message) {
  _messages.push_back(message);
  return true;
}

bool MessageBox::addText(std::string message) {
  return false;
}

void MessageBox::showText(bool whether_blank) {
  if (whether_blank) {
    for (auto & _message : _messages) {
      std::cout << _message << std::endl;
      std::cout << "====================================" << std::endl;
    }
  } else {
    for (auto & _message : _messages) {
      std::cout << _message << std::endl;
    }
  }
}




