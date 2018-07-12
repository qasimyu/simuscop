// ***************************************************************************
// TheadPool.h (c) 2015 Zhenhua Yu <yzh163@mail.ustc.edu.cn>
// Health Informatics Lab, University of Science and Technology of China
// All rights reserved.

#ifndef _THREADPOOL_H
#define _THREADPOOL_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
#include <cerrno>
#include <cassert>
#include <cstring>
#include <chrono>
#include <random>
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>

using namespace std;


class Thread {
	private:
		pthread_t tid;

		enum ThreadState {
			FREE = 0,
			BUSY = 1
		};
		mutable pthread_mutex_t pm;
		ThreadState state;
	public:
		Thread() {
			state = FREE;
			tid = pthread_self();
			pthread_mutex_init(&pm, NULL);
		}
		
		bool IsBusy() {return state == BUSY;}
		
		void Start() {
			pthread_mutex_lock(&pm);
			state = BUSY;
			pthread_mutex_unlock(&pm);
		}
		
		void Finish() {
			pthread_mutex_lock(&pm);
			state = FREE;
			pthread_mutex_unlock(&pm);
		}
		
		bool waitFinished() {
			usleep(50000);
			pthread_mutex_lock(&pm);
			if(state == FREE) {
				pthread_mutex_unlock(&pm);
				return true;
			}
			pthread_mutex_unlock(&pm);
			return false;
		}
		
		pthread_t getThreadId() const {return tid;}
		pthread_t* getThreadEntrance() {return &tid;}
};

class Work {
	private:
		void* (*process) (void* arg);
		void* arg;
		int wid;
		pthread_t tid;
		
		enum WorkState {
			READY = 0,
			UNDER_PROCESS = 1,
			FINISHED = 2
		};
		
		WorkState state;
		mutable pthread_mutex_t pm;
		struct timeval tv_started;
		struct timeval tv_finished;
	public:
		Work(void* (*process) (void* arg), void* arg, int wid) : process(process), arg(arg), wid(wid) {
			state = READY;
			//tid = 0;
			memset(&tv_started, 0, sizeof(tv_started));
			memset(&tv_finished, 0, sizeof(tv_finished));
			pthread_mutex_init(&pm, NULL);
		}
		
		bool IsReady() {
			return state == READY;
		}
		
		bool IsUnderProcess() {
			return state == UNDER_PROCESS;
		}
		
		bool IsFinished() {
			return state == FINISHED;
		}
		
		void StartProcess(pthread_t tid) {
			pthread_mutex_lock(&pm);
			this->tid = tid;
			gettimeofday(&tv_started, NULL);
			state = UNDER_PROCESS;
			pthread_mutex_unlock(&pm);
		}
		
		void FinishProcess() {
			pthread_mutex_lock(&pm);
			gettimeofday(&tv_finished, NULL);
			state = FINISHED;
			pthread_mutex_unlock(&pm);
		}
		
		long long ElapsedTime() const {
			pthread_mutex_lock(&pm);
			long long tv_usec;
			if(state == FINISHED) {
				tv_usec = (unsigned long long) (tv_finished.tv_sec-tv_started.tv_sec)*1000000 + (tv_finished.tv_usec - tv_started.tv_usec);
			}
			else {
				tv_usec = -1;
			}
			pthread_mutex_unlock(&pm);
			return tv_usec;
		}
		
		int getWorkId() const {return wid;}
		
		void run() {(*process)(arg);}
};

typedef struct queue {
	Work *work;
	struct queue *next;
}work_queue;

class ThreadPool {
	private:
		int max_thread_num;
		
		vector<Thread*> thread_list;
		work_queue *work_list;
		work_queue *rear;
		int cur_queue_size;
		int finishedWorks;
		bool shutdown;
		
		pthread_mutex_t work_lock;
		pthread_cond_t work_ready;
		
		long minRandNumber, maxRandNumber;
		map<pthread_t, default_random_engine> realGenerators;
		map<pthread_t, default_random_engine> intGenerators;
		
		void thread_routine();
		static void* threadFun(void *object) {
			((ThreadPool *)object)->thread_routine();
			return NULL;
		}
		
	public:
		ThreadPool(int max_thread_num) : max_thread_num(max_thread_num) {}
		
		ThreadPool() {	
			max_thread_num = 1;
		}
		
		void setThreadNumber(int threads) {
			max_thread_num = threads;
		}
		
		~ThreadPool();
		
		int worksFinished() {return finishedWorks;}
		
		void pool_init();	
		
		void pool_destroy();
		
		void clearWorks();
		
		void pool_add_work(void *(*process) (void *arg), void *arg, int wid);
		
		void wait();
		
		double randomDouble(double start, double end);
		long randomInteger(long start, long end);
		
};

#endif

