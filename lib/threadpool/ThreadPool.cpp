// ***************************************************************************
// TheadPool.cpp (c) 2015 Zhenhua Yu <yzh163@mail.ustc.edu.cn>
// Health Informatics Lab, University of Science and Technology of China
// All rights reserved.

#include "ThreadPool.h"

void ThreadPool::pool_init() {
	pthread_mutex_init(&work_lock, NULL); 
	pthread_cond_init(&work_ready, NULL);
	cur_queue_size = 0;
	finishedWorks = 0;
	shutdown = false;
	work_list = (work_queue*) malloc(sizeof(work_queue));
	work_list->next = NULL;
	rear = work_list;

	int thread_num = sysconf(_SC_NPROCESSORS_CONF);
	//cerr << "number of threads: " << thread_num << endl;

	max_thread_num = min(max_thread_num, thread_num);

	cerr << "\nnumber of available threads: " << thread_num << endl;
	cerr << "number of used threads: " << max_thread_num << endl << endl;

	for(size_t i = 0; i < max_thread_num; i++) {
		Thread* thread = new Thread();
		
		#ifdef __linux__
			cpu_set_t cpuset;
			CPU_ZERO(&cpuset);
			CPU_SET(i, &cpuset);
			if (pthread_setaffinity_np(thread->getThreadId(), sizeof(cpuset), &cpuset) < 0) {
				cerr << "set thread affinity failed." << endl;
			}
		#endif
		
		pthread_create(thread->getThreadEntrance(), NULL, threadFun, this);
		thread_list.push_back(thread);
		
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator(seed);
		minRandNumber = generator.min();
		maxRandNumber = generator.max();
		//uniform_real_distribution<double> realDist(0, maxNumber);
		//uniform_int_distribution<long> intDist(0, (long) maxNumber);
		realGenerators.insert(make_pair(thread->getThreadId(), generator));
		intGenerators.insert(make_pair(thread->getThreadId(), generator));
	}
}

void ThreadPool::pool_destroy() {
	if(thread_list.size() == 0 && work_list == NULL) {
		return;
	}
	shutdown = true;
	pthread_cond_broadcast(&work_ready);
	
	for(size_t i = 0; i < max_thread_num; i++) {
		pthread_join(thread_list[i]->getThreadId(), NULL);
	}
	
	std::vector<Thread*>::iterator it, it_e = thread_list.end();
	for(it = thread_list.begin(); it != it_e; it++) {
		delete *it;
	}
	thread_list.clear();
	
	work_queue *q,*p = work_list->next;
	while(p) {
		q = p->next;
		delete p->work;
		free(p);
		p = q;
	}
	free(work_list);
	work_list = NULL;
	rear = NULL;
	
	pthread_mutex_destroy(&work_lock);  
	pthread_cond_destroy(&work_ready);
	
}

void ThreadPool::pool_add_work(void *(*process) (void *arg), void *arg, int wid) {					
	Work *newwork = new Work(process, arg, wid);
	
	work_queue *p = (work_queue*) malloc(sizeof(work_queue));
	p->work = newwork;
	p->next = NULL;
	
	pthread_mutex_lock(&work_lock);
	rear->next = p;
	rear = p;
	
	assert(work_list->next != NULL);
	
	cur_queue_size++;
	
	pthread_mutex_unlock(&work_lock);
	pthread_cond_signal(&work_ready);
}

void ThreadPool::thread_routine() {
	//printf("starting thread 0x%x\n", pthread_self());
	while(1) {
		pthread_mutex_lock(&work_lock);
		while(cur_queue_size == 0) {
			//printf("thread 0x%x is waiting\n", pthread_self());
			if(shutdown) {
				pthread_mutex_unlock(&work_lock);
				pthread_exit(NULL);
			}
			pthread_cond_wait(&work_ready, &work_lock);
		}
		
		assert(cur_queue_size != 0);
		assert(work_list->next != NULL);
		
		size_t i, j;
		for(i = 0; i < thread_list.size(); i++) {
			if(thread_list[i]->getThreadId() == pthread_self()) {
				thread_list[i]->Start();
				break;
			}
		}
		
		work_queue *q = work_list, *p;
		p = q->next;
		while(p) {
			if(p->work->IsReady()) {
				break;
			}
			q = p;
			p = p->next;
		}
		
		assert(p != NULL);
		
		cur_queue_size--;
		q->next = p->next;
		if(rear == p) {
			rear = q;
		}
		
		pthread_mutex_unlock(&work_lock);
		
		//printf("thread 0x%x is starting to work\n", pthread_self());
		p->work->StartProcess(pthread_self());
		p->work->run();
		
		thread_list[i]->Finish();
		p->work->FinishProcess();
		long long tv_usec = p->work->ElapsedTime();
		int wid = p->work->getWorkId();
		
		delete p->work;
		free(p);
		p = NULL;
		
		finishedWorks++;
		
		
		//printf("work %d is finished in %d seconds by thread 0x%x\n", wid, tv_usec/1000000, pthread_self());
		
	}
}

void ThreadPool::clearWorks() {
	pthread_mutex_lock(&work_lock);
	
	work_queue *q, *p = work_list->next;
	while(p) {
		q = p->next;
		delete p->work;
		free(p);
		p = q;
	}
	work_list->next = NULL;
	cur_queue_size = 0;
	finishedWorks = 0;
	pthread_mutex_unlock(&work_lock);
}

void ThreadPool::wait() {
	while(1) {
		if(cur_queue_size == 0) {
			std::vector<Thread*>::iterator it;
			for(it = thread_list.begin(); it != thread_list.end(); it++) {
				if((*it)->IsBusy()) {
					break;
				}
			}
			if(it == thread_list.end()) {
				break;
			}
		}
		usleep(50000);
	}
}

double ThreadPool::randomDouble(double start, double end) {
	pthread_t tid = pthread_self();
	double number = realGenerators[tid]();
	return start+(end-start)*((number-minRandNumber)/(maxRandNumber-minRandNumber+1.0));
}
long ThreadPool::randomInteger(long start, long end) {
	pthread_t tid = pthread_self();
	double number = intGenerators[tid]();
	/*
	char cmd[100];
	sprintf(cmd, "echo %ld >> %x.int", number, tid);
	system(cmd);
	*/
	return start+(end-start)*((number-minRandNumber)/(maxRandNumber-minRandNumber+1.0));
}

ThreadPool::~ThreadPool() {
	std::vector<Thread*>::iterator it, it_e = thread_list.end();
	for(it = thread_list.begin(); it != it_e; it++) {
		delete *it;
	}
	thread_list.clear();

	work_queue *p = work_list->next;
	while(work_list) {
		delete work_list->work;
		free(work_list);
		work_list = p;
		p = work_list->next;
	}

	pthread_mutex_destroy(&work_lock);  
	pthread_cond_destroy(&work_ready);
}



