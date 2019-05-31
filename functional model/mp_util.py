from multiprocessing import Process, Semaphore, Queue

def run_process(semaphore, result_queue, index, func, input_value):
    result = func(input_value)
    # Release semaphore before queue put to allow main process to receive from
    #queue while we're writing.
    semaphore.release()
    result_queue.put((index, result))


class UniqueProcessMap(object):
    """A multi-processing version of map which uses a new process for each job,
    rather than worker processes that are re-used."""

    def __init__(self, max_processes):
        self.semaphore = Semaphore(0)
        self.result_queue = Queue()
        self.max_processes = max_processes

    def map(self, func, inputs):
        inputs = list(inputs)
        num_unfinished = len(inputs)
        results = [None] * len(inputs)

        for i in range(min(self.max_processes, len(inputs))):
            input_index = len(inputs) - 1
            process = Process(target=run_process, 
                              args=(self.semaphore, 
                                    self.result_queue, 
                                    input_index, func, 
                                    inputs.pop()))            
            process.start()

        while num_unfinished > 0:
            self.semaphore.acquire()
            result_index, result = self.result_queue.get()
            results[result_index] = result
            num_unfinished -= 1
            if len(inputs) > 0:
                input_index = len(inputs) - 1
                process = Process(target=run_process, 
                                  args=(self.semaphore, 
                                        self.result_queue, 
                                        input_index, func, 
                                        inputs.pop()))            
                process.start()

        return results
            


