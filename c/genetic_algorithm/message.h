//
// Created by 32270 on 2022/3/8.
//

#ifndef GENETIC_ALGORITHM_MESSAGE_H
#define GENETIC_ALGORITHM_MESSAGE_H

#include <vector>
#include <iostream>

class MessageBox{
private:
  std::vector<char *> _messages;
public:
  bool addText(char *message);
  bool addText(std::string message);

  void showText(bool whether_blank = false);
};

#endif //GENETIC_ALGORITHM_MESSAGE_H
