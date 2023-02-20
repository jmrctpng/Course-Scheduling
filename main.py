from typing import List, Tuple, Any
import random


class Data:
    professor = []
    info_technology = []
    comp_science = []
    lab_room = []
    reg_room = []
    lab_course = []
    reg_course = []

    @staticmethod
    def add_lab_course(course_code: str, lab_hours: int, lecture_hours: int, available_rooms: List[str]):
        course = [course_code, [lab_hours, lecture_hours], available_rooms]
        Data.lab_course.append(course)

    @staticmethod
    def add_lab_room(code: str, capacity: int):
        Data.lab_room.append((code, capacity))

    @staticmethod
    def add_reg_room(code: str, capacity: int):
        Data.reg_room.append((code, capacity))


class GeneticAlgorithm:
    chromosome = []

    @staticmethod
    def encode_chromosome() -> list[Any]:
        """
        Randomly select one item from each population list and create a new chromosome
        list containing those four items. The for loop iterates over the length of
        self.course_population and creates one chromosome for each course.
        """

        population_class = Population()

        for i in range(len(population_class.courses)):
            class_slot = population_class.class_slots.pop(random.randrange(len(population_class.class_slots)))
            course = population_class.courses.pop(random.randrange(len(population_class.courses)))
            block = population_class.blocks.pop(random.randrange(len(population_class.blocks)))
            professor = population_class.professors.pop(random.randrange(len(population_class.professors)))

            gene = Gene(class_slot, course, block, professor)

            GeneticAlgorithm.chromosome.append(gene)

    @staticmethod
    def calculate_fit():
        pass


class Gene:
    def __init__(self, class_slot, course, block, professor):
        self.class_sot = class_slot
        self.course = course
        self.block = block
        self.professor = professor

    def get_class_slot(self):
        return self.class_sot

    def get_course(self):
        return self.course

    def get_block(self):
        return self.block

    def get_professor(self):
        return self.professor()

    def get_fitness_score(self):
        pass


class GeneFitness:

    def __init__(self, gene):
        self.gene = gene

    def classroom_capacity(self):
        if self.gene.block().get_enrolled_students() <= self.gene.class_slot()[0].get_capacity():
            return 1
        else:
            return 0

    def professor_work_load(self):
        day = self.gene.get_class_slot()[1]
        day_sched = []
        schedule = self.gene.professor.get_schedule()

        if day == 'M':
            day_sched = schedule.get_monday()
        elif day == 'T':
            day_sched = schedule.get_tuesday()
        elif day == 'W':
            day_sched = schedule.get_wednesday()
        elif day == 'Th':
            day_sched = schedule.get_thursday()
        elif day == 'F':
            day_sched = schedule.get_friday()
        elif day == 'S':
            day_sched = schedule.get_saturday()

        result = has_overlap(self.gene.class_slot[2][0], day_sched)

        return -1 if result else 0

    def schedule_availability(self):
        pass


class ClassRoom:

    def __init__(self, code, capacity, schedule):
        self.code = code
        self.capacity = capacity
        self.schedule = schedule

    def get_code(self):
        return self.code

    def get_capacity(self):
        return self.capacity()

    def get_schedule(self):
        return self.schedule

class Course:
    lab_courses = []
    reg_courses = []

    def __init__(self, code, lab_hours, lecture_hours, available_rooms):
        self.code = code
        self.lab_hours = lab_hours
        self.lecture_hours = lecture_hours
        self.available_rooms = available_rooms

    def get_code(self):
        return self.code()

    def get_lab_hours(self):
        return self.lab_hours()

    def get_lecture_hours(self):
        return self.lecture_hours()

    # --------------------------------------------------------------------------------------------------

    @staticmethod
    def add_lab_course(course_code: str, lab_hours: int, lecture_hours: int, available_rooms: List[str]):
        course = [course_code, [lab_hours, lecture_hours], available_rooms]
        Course.lab_courses.append(course)

    @staticmethod
    def add_reg_course(course_code: str, lab_hours: int, lecture_hours: int, available_rooms: List[str]):
        course = [course_code, [lab_hours, lecture_hours], available_rooms]
        Course.reg_courses.append(course)


class Schedule:
    def __init__(self):
        self.monday = []
        self.tuesday = []
        self.wednesday = []
        self.thursday = []
        self.friday = []
        self.saturday = []

    def get_monday(self):
        return self.monday

    def set_monday(self, monday):
        self.monday = monday

    def get_tuesday(self):
        return self.tuesday

    def set_tuesday(self, tuesday):
        self.tuesday = tuesday

    def get_wednesday(self):
        return self.wednesday

    def set_wednesday(self, wednesday):
        self.wednesday = wednesday

    def get_thursday(self):
        return self.thursday

    def set_thursday(self, thursday):
        self.thursday = thursday

    def get_friday(self):
        return self.friday

    def set_friday(self, friday):
        self.friday = friday

    def get_saturday(self):
        return self.saturday

    def set_saturday(self, saturday):
        self.saturday = saturday


class Professor:
    professors = []

    def __init__(self, prof_id, course_code_handle, schedule):
        self.prof_id = prof_id
        self.course_code_handle = course_code_handle
        self.schedule = schedule

    def get_prof_id(self):
        return self.prof_id()

    def get_course_code_handle(self):
        return self.course_code_handle()

    def get_schedule(self):
        return self.schedule()

    # -------------------------------------------------
    @staticmethod
    def add_prof(prof_id: str, course_code_handle: str):
        element = [prof_id, course_code_handle]
        Professor.professors.append(element)


class Program:

    def __init__(self, block_code, enrolled_students, courses):
        self.block_code = block_code
        self.enrolled_students = enrolled_students
        self.courses = courses

    def get_block_code(self):
        return self.block_code

    def get_enrolled_students(self):
        return self.enrolled_students

    def get_courses(self):
        return self.courses


class Population:
    def __init__(self):
        self.population = []
        self.courses = []
        self.blocks = []
        self.professors = []
        self.class_slots = []

    def initialize_population(self, classSlot, professor, course, block):
        pass

    def generate_course_population(self):
        block = Data()
        available_rooms = block.comp_science + block.info_technology

        for room in available_rooms:
            for course in room[3]:
                self.courses.append(course)

    def generate_block_population(self):
        block = Data()
        available_rooms = block.comp_science + block.info_technology

        for room in available_rooms:
            for i in range(len(room[3])):
                self.blocks.append(room[i][0])

    def generate_prof_population(self):
        professor = Data()

        for prof in professor.professor:
            for i in range(len(prof[1])):
                self.professors.append(prof[0])

    def generate_class_slot_population(self, meetingTime, day):
        class_slot = labRoomMapping(meetingTime, day) + regRoomMapping(day)
        self.class_slots.append(class_slot)

    def get_population(self):
        return self.population


def labRoomMapping(meetingTime: List[str], meetingDays: List[str]) -> List[List[str]]:
    """
        This function maps the available laboratory rooms with the given meeting time and meeting days.

        Arguments:
        meetingTime -- A list of meeting hours in the format "HH".
        meetingDays -- A list of meeting days.

        Returns:
        lab_map -- A list of lists, where each inner list represents the mapping of a laboratory room with a day and three hours of meeting time.

        Time complexity:  O(len(labRoomList) * (len(meetingTime)/3) * len(meetingDays)) or O(n^3)
        This is because, for each laboratory room in labRoomList, it loops through each day in meetingDays and each three-hour time slot in meetingTime. The length of each time slot is determined by the range function, where it steps by 3. The time complexity is proportional to the number of iterations required to map each laboratory room with its corresponding day and time slot.
        """

    labRoom = Data()
    lab_map = [[lab, day, meetingTime[i:i + 4]]
               for lab in labRoom.lab_room
               for day in meetingDays
               for i in range(0, len(meetingTime), 3)
               if len(meetingTime[i:i + 4]) == 4]

    return lab_map


def regRoomMapping(meetingDays: List[str]) -> List[List[str]]:
    reg_room = Data()
    pattern = [["07", "08", "09"], ["07", "08", "09", "10"], ["09", "10", "11"], ["09", "10", "11", "12"],
               ["10", "11", "12"], ["10", "11", "12", "13"],
               ["11", "12", "13"], ["11", "12", "13", "14"], ["12", "13", "14"], ["12", "13", "14", "15"],
               ["13", "14", "15"], ["13", "14", "15", "16"],
               ["14", "15", "16"], ["14", "15", "16", "17"], ["15", "16", "17"], ["16", "17", "18", "19"],
               ["17", "18", "19"]]

    reg_room_map = [[room, day, meetingTime]
                    for room in reg_room.reg_room
                    for day in meetingDays
                    for meetingTime in pattern
                    ]

    return reg_room_map


def has_overlap(time_slot, schedule):
    for slot in time_slot:
        for i, time in enumerate(schedule):
            if slot[0] in time and i < len(schedule) - len(slot) + 1:
                if all(slot[j] == schedule[i+j] for j in range(1, len(slot))):
                    return True
    return False


if __name__ == '__main__':
    meetingTime = ["07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
    day = ["M", "T", "W", "Th", "F", "S"]

    Data.add_lab_room("SL1", 50)
    Data.add_lab_room("SL2", 50)
    Data.add_lab_room("SL3", 50)
    Data.add_lab_room("CISCO", 50)
    Data.add_lab_room("ML", 50)
    Data.add_lab_room("ITL", 50)
    Data.add_reg_room("101", 50)
    Data.add_reg_room("102", 50)
    Data.add_reg_room("103", 50)
    Data.add_reg_room("104", 50)
    Data.add_reg_room("201", 50)
    Data.add_reg_room("202", 50)
    Data.add_reg_room("203", 50)
    Data.add_reg_room("EDL", 50)
    Data.add_reg_room("401", 50)
    Data.add_reg_room("402", 50)
    Data.add_reg_room("403", 50)

    Data.add_lab_course("CS 121", 2, 3, ["Sl1", "SL2", "SL3"])
    Data.add_lab_course("IT 221", 2, 3,
                          ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                           "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("CS 222", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    Data.add_lab_course("CS 221", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    Data.add_lab_course("CS 322", 2, 3, ["SL1", "SL2", "SL3"])
    Data.add_lab_course("CS 324", 2, 3,
                          ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                           "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("CSE 401", 2, 3,
                          ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                           "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("CS 422", 2, 3, ["ML"])
    Data.add_lab_course("IT 323", 2, 3, ["ML"])
    Data.add_lab_course("CSE 403", 2, 3, ["ML"])
    Data.add_lab_course("CS 111", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    Data.add_lab_course("IT 221", 2, 3,
                          ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                           "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("IT 223", 2, 3, ["CISCO"])
    Data.add_lab_course("IT 222", 2, 3,
                          ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                           "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("NTT 403", 2, 3, ["CISCO"])
    Data.add_lab_course("NTT 404", 2, 3, ["CISCO"])
    Data.add_lab_course("IT 322", 2, 3,
                          ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                           "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("IT 323", 2, 3,
                          ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                           "ML", "ITL", "EDL", "401", "402", "403"])
