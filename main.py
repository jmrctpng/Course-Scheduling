from typing import List, Tuple


class ClassRoom:
    def __init__(self):
        self.lab_rooms = []
        self.reg_rooms = []

    def add_lab_room(self, code: str, capacity: int):
        self.lab_rooms.append((code, capacity))

    def add_reg_room(self, code: str, capacity: int):
        self.reg_rooms.append((code, capacity))

    def show_all_lab_rooms(self) -> List[Tuple[str, int]]:
        return self.lab_rooms

    def show_all_reg_rooms(self) -> List[Tuple[str, int]]:
        return self.reg_rooms


class Course:
    lab_courses = []
    reg_courses = []

    @staticmethod
    def add_lab_course(course_code: str, lab_hours: int, lecture_hours: int, available_rooms: List[str]):
        course = [course_code, [lab_hours, lecture_hours], available_rooms]
        Course.lab_courses.append(course)

    @staticmethod
    def add_reg_course(course_code: str, lab_hours: int, lecture_hours: int, available_rooms: List[str]):
        course = [course_code, [lab_hours, lecture_hours], available_rooms]
        Course.reg_courses.append(course)


class Professor:
    def __init__(self):
        self.prof_list = []

    def add_prof(self, prof_id: str, course_code_handle: str):
        element = [prof_id, course_code_handle]
        self.prof_list.append(element)


class Program:
    def __init__(self):
        self.info_technology_list = []
        self.comp_sci_list = []

    def add_info_tech_block(self, block_code: str, total_students: int):
        """
        Adds a block of information technology to the program.
        """
        self.info_technology_list.append((block_code, total_students))

    def add_comp_sci_block(self, block_code: str, total_students: int):
        """
        Adds a block of computer science to the program.
        """
        self.comp_sci_list.append((block_code, total_students))


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

    labRoom = ClassRoom()
    lab_map = [[lab, day, meetingTime[i:i + 4]]
               for lab in labRoom.lab_rooms
               for day in meetingDays
               for i in range(0, len(meetingTime), 3)
               if len(meetingTime[i:i + 4]) == 4]

    return lab_map


def regRoomMapping(meetingTime: List[str], meetingDays: List[str]) -> List[List[str]]:
    reg_room_map = []
    reg_room = ClassRoom()

    #not complete

    return reg_room_map


if __name__ == '__main__':

    meetingTime = ["07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
    day = ["M", "T", "W", "Th", "F", "S"]

    available_room = ClassRoom()
    available_room.add_lab_room("SL1", 50)
    available_room.add_lab_room("SL2", 50)
    available_room.add_lab_room("SL3", 50)
    available_room.add_lab_room("CISCO", 50)
    available_room.add_lab_room("ML", 50)
    available_room.add_lab_room("ITL", 50)
    available_room.add_reg_room("101", 50)
    available_room.add_reg_room("102", 50)
    available_room.add_reg_room("103", 50)
    available_room.add_reg_room("104", 50)
    available_room.add_reg_room("201", 50)
    available_room.add_reg_room("202", 50)
    available_room.add_reg_room("203", 50)
    available_room.add_reg_room("EDL", 50)
    available_room.add_reg_room("401", 50)
    available_room.add_reg_room("402", 50)
    available_room.add_reg_room("403", 50)

    available_course = Course()
    available_course.add_lab_course("CS 121", 2, 3, ["Sl1", "SL2", "SL3"])
    available_course.add_lab_course("IT 221", 2, 3, ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO", "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("CS 222", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    available_course.add_lab_course("CS 221", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    available_course.add_lab_course("CS 322", 2, 3, ["SL1", "SL2", "SL3"])
    available_course.add_lab_course("CS 324", 2, 3, ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO", "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("CSE 401", 2, 3, ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO", "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("CS 422", 2, 3, ["ML"])
    available_course.add_lab_course("IT 323", 2, 3, ["ML"])
    available_course.add_lab_course("CSE 403", 2, 3, ["ML"])
    available_course.add_lab_course("CS 111", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    available_course.add_lab_course("IT 221", 2, 3, ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO", "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("IT 223", 2, 3, ["CISCO"])
    available_course.add_lab_course("IT 222", 2, 3, ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO", "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("NTT 403", 2, 3, ["CISCO"])
    available_course.add_lab_course("NTT 404", 2, 3, ["CISCO"])
    available_course.add_lab_course("IT 322", 2, 3, ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO", "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("IT 323", 2, 3, ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO", "ML", "ITL", "EDL", "401", "402", "403"])

