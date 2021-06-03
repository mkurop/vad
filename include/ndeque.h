//
// Created by marcin on 11.01.19.
//

#ifndef ZOOMER_NDEQUE_H
#define ZOOMER_NDEQUE_H

#include <deque>
#include <iostream>
template<class DequeElement>
class NDeque {
  std::deque<DequeElement> store;
  size_t maxSize{0};
 public:
  NDeque() : maxSize(0) {}
  NDeque(size_t maxSize) : maxSize(maxSize) {}

  void push_back(const DequeElement &dequeElement) {
    store.push_back(dequeElement);
    if (store.size() > maxSize) {
      store.pop_front();
    }
  }

  void push_front(const DequeElement &dequeElement) {
    store.push_front(dequeElement);
    if (store.size() > maxSize) {
      store.pop_back();
    }
  }

  DequeElement &back() {
    return store.back();
  }

  DequeElement &front() {
    return store.front();
  }

  template<class... Args>
  void emplace_back(Args &&... args) {
    store.emplace_back(std::forward<Args>(args)...);
    if (store.size() > maxSize) {
      store.pop_front();
    }
  }

  template<class... Args>
  void emplace_front(Args &&... args) {
    store.emplace_front(std::forward<Args>(args)...);
    if (store.size() > maxSize) {
      store.pop_back();
    }
  }

  void resize(size_t newMaxSize) {
    maxSize = newMaxSize;
  }

  auto begin() -> decltype(store.begin()) {
    return store.begin();
  }

  auto end() -> decltype(store.end()) {
    return store.end();
  }

  auto rend() -> decltype(store.rend()) {
      return store.rend();
  }

  auto rbegin() -> decltype(store.rbegin()) {
      return store.rbegin();
  }

  bool empty() const {
    return store.empty();
  }

  void clear() {
    store.clear();
  }

  void pop_front() {
    store.pop_front();
  }

  void pop_back() {
    store.pop_back();
  }

  size_t size() const {
    return store.size();
  }
  size_t size() {
    return store.size();
  }
};

#endif //ZOOMER_NDEQUE_H
