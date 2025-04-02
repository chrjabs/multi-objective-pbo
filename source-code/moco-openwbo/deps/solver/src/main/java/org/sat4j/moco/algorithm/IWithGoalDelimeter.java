package org.sat4j.moco.algorithm;

import org.sat4j.moco.goal_delimeter.GoalDelimeterI;


interface IWithGoalDelimeter {
    
    void setGoalDelimeter(GoalDelimeterI gd);
    GoalDelimeterI GetGoalDelimeter();
}
