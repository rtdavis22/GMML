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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "gmml/internal/stubs/logging.h"

#include <stdio.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN  // We only need minimal includes
#include <windows.h>
#define snprintf _snprintf    // see comment in strutil.cc
#elif defined(HAVE_PTHREAD)
#include <pthread.h>
#endif

using std::string;

namespace gmml {
namespace internal {

// ===================================================================
// emulates google3/base/mutex.cc

#ifdef _WIN32

struct Mutex::Internal {
  CRITICAL_SECTION mutex;
#ifndef NDEBUG
  // Used only to implement AssertHeld().
  DWORD thread_id;
#endif
};

Mutex::Mutex()
  : mInternal(new Internal) {
  InitializeCriticalSection(&mInternal->mutex);
}

Mutex::~Mutex() {
  DeleteCriticalSection(&mInternal->mutex);
  delete mInternal;
}

void Mutex::Lock() {
  EnterCriticalSection(&mInternal->mutex);
#ifndef NDEBUG
  mInternal->thread_id = GetCurrentThreadId();
#endif
}

void Mutex::Unlock() {
#ifndef NDEBUG
  mInternal->thread_id = 0;
#endif
  LeaveCriticalSection(&mInternal->mutex);
}

void Mutex::AssertHeld() {
#ifndef NDEBUG
  GOOGLE_DCHECK_EQ(mInternal->thread_id, GetCurrentThreadId());
#endif
}

#elif defined(HAVE_PTHREAD)

struct Mutex::Internal {
  pthread_mutex_t mutex;
};

Mutex::Mutex()
  : mInternal(new Internal) {
  pthread_mutex_init(&mInternal->mutex, NULL);
}

Mutex::~Mutex() {
  pthread_mutex_destroy(&mInternal->mutex);
  delete mInternal;
}

void Mutex::Lock() {
  int result = pthread_mutex_lock(&mInternal->mutex);
  if (result != 0) {
    GOOGLE_LOG(FATAL) << "pthread_mutex_lock: " << strerror(result);
  }
}

void Mutex::Unlock() {
  int result = pthread_mutex_unlock(&mInternal->mutex);
  if (result != 0) {
    GOOGLE_LOG(FATAL) << "pthread_mutex_unlock: " << strerror(result);
  }
}

void Mutex::AssertHeld() {
  // pthreads dosn't provide a way to check which thread holds the mutex.
  // TODO(kenton):  Maybe keep track of locking thread ID like with WIN32?
}

#endif

void DefaultLogHandler(LogLevel level, const char* filename, int line,
                       const string& message) {
  static const char* level_names[] = { "INFO", "WARNING", "ERROR", "FATAL" };

  // We use fprintf() instead of cerr because we want this to work at static
  // initialization time.
  fprintf(stderr, "[GMML %s %s:%d] %s\n",
          level_names[level], filename, line, message.c_str());
  fflush(stderr);  // Needed on MSVC.
}

void NullLogHandler(LogLevel level, const char* filename, int line,
                    const string& message) {
  // Nothing.
}

// I may want to either keep this as the default log handler, use
// DefaultLogHandler above, or write a log handler that writes to a file
// depending on how logging is used in the library.
static LogHandler* log_handler_ = &DefaultLogHandler;

/*
static int log_silencer_count_ = 0;

static Mutex* log_silencer_count_mutex_ = NULL;
GOOGLE_PROTOBUF_DECLARE_ONCE(log_silencer_count_init_);

void DeleteLogSilencerCount() {
  delete log_silencer_count_mutex_;
  log_silencer_count_mutex_ = NULL;
}
void InitLogSilencerCount() {
  log_silencer_count_mutex_ = new Mutex;
  //OnShutdown(&DeleteLogSilencerCount);
}
void InitLogSilencerCountOnce() {
  GoogleOnceInit(&log_silencer_count_init_, &InitLogSilencerCount);
}
*/

LogMessage::LogMessage(LogLevel level, const char* filename, int line)
  : level_(level), filename_(filename), line_(line) {}
LogMessage::~LogMessage() {}

LogMessage& LogMessage::operator<<(const string& value) {
  message_ += value;
  return *this;
}

LogMessage& LogMessage::operator<<(const char* value) {
  message_ += value;
  return *this;
}

// Since this is just for logging, we don't care if the current locale changes
// the results -- in fact, we probably prefer that.  So we use snprintf()
// instead of Simple*toa().
#undef DECLARE_STREAM_OPERATOR
#define DECLARE_STREAM_OPERATOR(TYPE, FORMAT)                       \
  LogMessage& LogMessage::operator<<(TYPE value) {                  \
    /* 128 bytes should be big enough for any of the primitive */   \
    /* values which we print with this, but well use snprintf() */  \
    /* anyway to be extra safe. */                                  \
    char buffer[128];                                               \
    snprintf(buffer, sizeof(buffer), FORMAT, value);                \
    /* Guard against broken MSVC snprintf(). */                     \
    buffer[sizeof(buffer)-1] = '\0';                                \
    message_ += buffer;                                             \
    return *this;                                                   \
  }

DECLARE_STREAM_OPERATOR(char         , "%c" )
DECLARE_STREAM_OPERATOR(int          , "%d" )
DECLARE_STREAM_OPERATOR(unsigned int , "%u" )
DECLARE_STREAM_OPERATOR(long         , "%ld")
DECLARE_STREAM_OPERATOR(unsigned long, "%lu")
DECLARE_STREAM_OPERATOR(double       , "%g" )
#undef DECLARE_STREAM_OPERATOR

void LogMessage::Finish() {
  bool suppress = false;

  if (level_ != LOGLEVEL_FATAL) {
    //InitLogSilencerCountOnce();
    //MutexLock lock(log_silencer_count_mutex_);
    //suppress = internal::log_silencer_count_ > 0;
  }

  if (!suppress) {
    internal::log_handler_(level_, filename_, line_, message_);
  }

  if (level_ == LOGLEVEL_FATAL) {
#ifdef PROTOBUF_USE_EXCEPTIONS
    throw FatalException(filename_, line_, message_);
#else
    abort();
#endif
  }
}

void LogFinisher::operator=(LogMessage& other) {
  other.Finish();
}

}  // namespace internal

LogHandler* SetLogHandler(LogHandler* new_func) {
  LogHandler* old = internal::log_handler_;
  if (old == &internal::NullLogHandler) {
    old = NULL;
  }
  if (new_func == NULL) {
    internal::log_handler_ = &internal::NullLogHandler;
  } else {
    internal::log_handler_ = new_func;
  }
  return old;
}

/*
LogSilencer::LogSilencer() {
  internal::InitLogSilencerCountOnce();
  MutexLock lock(internal::log_silencer_count_mutex_);
  ++internal::log_silencer_count_;
};

LogSilencer::~LogSilencer() {
  internal::InitLogSilencerCountOnce();
  MutexLock lock(internal::log_silencer_count_mutex_);
  --internal::log_silencer_count_;
};



#ifdef _WIN32

struct ProtobufOnceInternal {
  ProtobufOnceInternal() {
    InitializeCriticalSection(&critical_section);
  }
  ~ProtobufOnceInternal() {
    DeleteCriticalSection(&critical_section);
  }
  CRITICAL_SECTION critical_section;
};

ProtobufOnceType::~ProtobufOnceType()
{
  delete internal_;
  internal_ = NULL;
}

ProtobufOnceType::ProtobufOnceType() {
  // internal_ may be non-NULL if Init() was already called.
  if (internal_ == NULL) internal_ = new ProtobufOnceInternal;
}

void ProtobufOnceType::Init(void (*init_func)()) {
  // internal_ may be NULL if we're still in dynamic initialization and the
  // constructor has not been called yet.  As mentioned in once.h, we assume
  // that the program is still single-threaded at this time, and therefore it
  // should be safe to initialize internal_ like so.
  if (internal_ == NULL) internal_ = new ProtobufOnceInternal;

  EnterCriticalSection(&internal_->critical_section);
  if (!initialized_) {
    init_func();
    initialized_ = true;
  }
  LeaveCriticalSection(&internal_->critical_section);
}

#endif
*/

}  // namespace gmml
