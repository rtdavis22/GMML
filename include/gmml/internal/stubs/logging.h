// Protocol Buffers - Google's data interchange format
// Copyright 2008 Google Inc.  All rights reserved.
// http://code.google.com/p/protobuf/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Note(Robert): I commented out the LogSilencer stuff because it uses pthreads
// and I don't want the user to have to link with a pthread library when there
// isn't any multithreading as of yet. The commented lines should eventually
// either be removed completely or uncommented.

#ifndef GMML_INTERNAL_STUBS_LOGGING_H_
#define GMML_INTERNAL_STUBS_LOGGING_H_

#include <stdlib.h>
#include <string.h>

#ifndef _WIN32
#include <pthread.h>
#endif

#include <string>

#include "gmml/internal/stubs/common.h"

namespace gmml {
namespace internal {

// A Mutex is a non-reentrant (aka non-recursive) mutex.  At most one thread T
// may hold a mutex at a given time.  If T attempts to Lock() the same Mutex
// while holding it, T will deadlock.
class Mutex {
 public:
  // Create a Mutex that is not held by anybody.
  Mutex();

  // Destructor
  ~Mutex();

  // Block if necessary until this Mutex is free, then acquire it exclusively.
  void Lock();

  // Release this Mutex.  Caller must hold it exclusively.
  void Unlock();

  // Crash if this Mutex is not held exclusively by this thread.
  // May fail to crash when it should; will never crash when it should not.
  void AssertHeld();

 private:
  struct Internal;
  Internal* mInternal;

  DISALLOW_COPY_AND_ASSIGN(Mutex);
};

// MutexLock(mu) acquires mu when constructed and releases it when destroyed.
class MutexLock {
 public:
  explicit MutexLock(Mutex *mu) : mu_(mu) { this->mu_->Lock(); }
  ~MutexLock() { this->mu_->Unlock(); }
 private:
  Mutex *const mu_;
  DISALLOW_COPY_AND_ASSIGN(MutexLock);
};

// TODO(kenton):  Implement these?  Hard to implement portably.
typedef MutexLock ReaderMutexLock;
typedef MutexLock WriterMutexLock;

// MutexLockMaybe is like MutexLock, but is a no-op when mu is NULL.
class MutexLockMaybe {
 public:
  explicit MutexLockMaybe(Mutex *mu) :
    mu_(mu) { if (this->mu_ != NULL) { this->mu_->Lock(); } }
  ~MutexLockMaybe() { if (this->mu_ != NULL) { this->mu_->Unlock(); } }
 private:
  Mutex *const mu_;
  DISALLOW_COPY_AND_ASSIGN(MutexLockMaybe);
};

}  // namespace internal

// We made these internal so that they would show up as such in the docs,
// but we don't want to stick "internal::" in front of them everywhere.
using internal::Mutex;
using internal::MutexLock;
using internal::ReaderMutexLock;
using internal::WriterMutexLock;
using internal::MutexLockMaybe;

// ===================================================================
// emulates google3/base/logging.h

enum LogLevel {
  LOGLEVEL_INFO,     // Informational.  This is never actually used by
                     // libprotobuf.
  LOGLEVEL_WARNING,  // Warns about issues that, although not technically a
                     // problem now, could cause problems in the future.  For
                     // example, a // warning will be printed when parsing a
                     // message that is near the message size limit.
  LOGLEVEL_ERROR,    // An error occurred which should never happen during
                     // normal use.
  LOGLEVEL_FATAL,    // An error occurred from which the library cannot
                     // recover.  This usually indicates a programming error
                     // in the code which calls the library, especially when
                     // compiled in debug mode.

#ifdef NDEBUG
  LOGLEVEL_DFATAL = LOGLEVEL_ERROR
#else
  LOGLEVEL_DFATAL = LOGLEVEL_FATAL
#endif
};

namespace internal {

class LogFinisher;

class LogMessage {
 public:
  LogMessage(LogLevel level, const char* filename, int line);
  ~LogMessage();

  LogMessage& operator<<(const std::string& value);
  LogMessage& operator<<(const char* value);
  LogMessage& operator<<(char value);
  LogMessage& operator<<(int value);
  LogMessage& operator<<(unsigned int value);
  LogMessage& operator<<(long value);
  LogMessage& operator<<(unsigned long value);
  LogMessage& operator<<(double value);

 private:
  friend class LogFinisher;
  void Finish();

  LogLevel level_;
  const char* filename_;
  int line_;
  std::string message_;
};

// Used to make the entire "LOG(BLAH) << etc." expression have a void return
// type and print a newline after each message.
class LogFinisher {
 public:
  void operator=(LogMessage& other);
};

}  // namespace internal

// Undef everything in case we're being mixed with some other Google library
// which already defined them itself.  Presumably all Google libraries will
// support the same syntax for these so it should not be a big deal if they
// end up using our definitions instead.
#undef GOOGLE_LOG
#undef GOOGLE_LOG_IF

#undef GOOGLE_CHECK
#undef GOOGLE_CHECK_EQ
#undef GOOGLE_CHECK_NE
#undef GOOGLE_CHECK_LT
#undef GOOGLE_CHECK_LE
#undef GOOGLE_CHECK_GT
#undef GOOGLE_CHECK_GE

#undef GOOGLE_DLOG
#undef GOOGLE_DCHECK
#undef GOOGLE_DCHECK_EQ
#undef GOOGLE_DCHECK_NE
#undef GOOGLE_DCHECK_LT
#undef GOOGLE_DCHECK_LE
#undef GOOGLE_DCHECK_GT
#undef GOOGLE_DCHECK_GE


#define GOOGLE_LOG(LEVEL)                                                 \
  ::gmml::internal::LogFinisher() =                           \
    ::gmml::internal::LogMessage(                             \
      ::gmml::LOGLEVEL_##LEVEL, __FILE__, __LINE__)
#define GOOGLE_LOG_IF(LEVEL, CONDITION) \
  !(CONDITION) ? (void)0 : GOOGLE_LOG(LEVEL)


#define GOOGLE_CHECK(EXPRESSION) \
  GOOGLE_LOG_IF(FATAL, !(EXPRESSION)) << "CHECK failed: " #EXPRESSION ": "
#define GOOGLE_CHECK_EQ(A, B) GOOGLE_CHECK((A) == (B))
#define GOOGLE_CHECK_NE(A, B) GOOGLE_CHECK((A) != (B))
#define GOOGLE_CHECK_LT(A, B) GOOGLE_CHECK((A) <  (B))
#define GOOGLE_CHECK_LE(A, B) GOOGLE_CHECK((A) <= (B))
#define GOOGLE_CHECK_GT(A, B) GOOGLE_CHECK((A) >  (B))
#define GOOGLE_CHECK_GE(A, B) GOOGLE_CHECK((A) >= (B))

#ifdef NDEBUG

#define GOOGLE_DLOG GOOGLE_LOG_IF(INFO, false)

#define GOOGLE_DCHECK(EXPRESSION) while(false) GOOGLE_CHECK(EXPRESSION)
#define GOOGLE_DCHECK_EQ(A, B) GOOGLE_DCHECK((A) == (B))
#define GOOGLE_DCHECK_NE(A, B) GOOGLE_DCHECK((A) != (B))
#define GOOGLE_DCHECK_LT(A, B) GOOGLE_DCHECK((A) <  (B))
#define GOOGLE_DCHECK_LE(A, B) GOOGLE_DCHECK((A) <= (B))
#define GOOGLE_DCHECK_GT(A, B) GOOGLE_DCHECK((A) >  (B))
#define GOOGLE_DCHECK_GE(A, B) GOOGLE_DCHECK((A) >= (B))

#else  // NDEBUG

#define GOOGLE_DLOG GOOGLE_LOG

#define GOOGLE_DCHECK    GOOGLE_CHECK
#define GOOGLE_DCHECK_EQ GOOGLE_CHECK_EQ
#define GOOGLE_DCHECK_NE GOOGLE_CHECK_NE
#define GOOGLE_DCHECK_LT GOOGLE_CHECK_LT
#define GOOGLE_DCHECK_LE GOOGLE_CHECK_LE
#define GOOGLE_DCHECK_GT GOOGLE_CHECK_GT
#define GOOGLE_DCHECK_GE GOOGLE_CHECK_GE

#endif  // !NDEBUG

typedef void LogHandler(LogLevel level, const char* filename, int line,
                        const std::string& message);

// The protobuf library sometimes writes warning and error messages to
// stderr.  These messages are primarily useful for developers, but may
// also help end users figure out a problem.  If you would prefer that
// these messages be sent somewhere other than stderr, call SetLogHandler()
// to set your own handler.  This returns the old handler.  Set the handler
// to NULL to ignore log messages (but see also LogSilencer, below).
//
// Obviously, SetLogHandler is not thread-safe.  You should only call it
// at initialization time, and probably not from library code.  If you
// simply want to suppress log messages temporarily (e.g. because you
// have some code that tends to trigger them frequently and you know
// the warnings are not important to you), use the LogSilencer class
// below.
LogHandler* SetLogHandler(LogHandler* new_func);

// Create a LogSilencer if you want to temporarily suppress all log
// messages.  As long as any LogSilencer objects exist, non-fatal
// log messages will be discarded (the current LogHandler will *not*
// be called).  Constructing a LogSilencer is thread-safe.  You may
// accidentally suppress log messages occurring in another thread, but
// since messages are generally for debugging purposes only, this isn't
// a big deal.  If you want to intercept log messages, use SetLogHandler().
/*
class LogSilencer {
 public:
  LogSilencer();
  ~LogSilencer();
};
*/
/*
#ifdef _WIN32

struct ProtobufOnceInternal;

struct ProtobufOnceType {
  ProtobufOnceType();
  ~ProtobufOnceType();
  void Init(void (*init_func)());

  volatile bool initialized_;
  ProtobufOnceInternal* internal_;
};

#define GOOGLE_PROTOBUF_DECLARE_ONCE(NAME)                    \
  ::gmml::ProtobufOnceType NAME

inline void GoogleOnceInit(ProtobufOnceType* once, void (*init_func)()) {
  // Note:  Double-checked locking is safe on x86.
  if (!once->initialized_) {
    once->Init(init_func);
  }
}

#else

typedef pthread_once_t ProtobufOnceType;

#define GOOGLE_PROTOBUF_DECLARE_ONCE(NAME)                    \
  pthread_once_t NAME = PTHREAD_ONCE_INIT

inline void GoogleOnceInit(ProtobufOnceType* once, void (*init_func)()) {
  pthread_once(once, init_func);
}

#endif
*/

}  // namespace gmml

#endif  // GMML_INTERNAL_STUBS_LOGGING_H_
