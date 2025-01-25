function [gen, det] = get_crc_objective(crc_length)
% 根据指定的 CRC 长度创建 CRC 生成器和检测器对象
%   crc_length - CRC 校验的长度（支持 4, 6, 8, 10, 12, 16, 24）
%   gen - CRC 生成器对象，用于为数据添加 CRC 校验码
%   det - CRC 检测器对象，用于检测数据的 CRC 校验码
    switch crc_length
        case 4
            gen = comm.CRCGenerator('Polynomial',[1 0 0 1 1],'InitialConditions',zeros(1, 4),'DirectMethod',true,'FinalXOR',zeros(1, 4));
            det = comm.CRCDetector('Polynomial',[1 0 0 1 1],'InitialConditions',zeros(1, 4),'DirectMethod',true,'FinalXOR',zeros(1, 4));
        case 6
            gen = comm.CRCGenerator('Polynomial',[1 0 0 0 0 1 1],'InitialConditions',zeros(1, 6),'DirectMethod',true,'FinalXOR',zeros(1, 6));
            det = comm.CRCDetector('Polynomial',[1 0 0 0 0 1 1],'InitialConditions',zeros(1, 6),'DirectMethod',true,'FinalXOR',zeros(1, 6));
        case 8
            gen = comm.CRCGenerator('Polynomial','0xA6','InitialConditions','0x00','DirectMethod',true,'FinalXOR','0x00');
            det = comm.CRCDetector('Polynomial','0xA6','InitialConditions','0x00','DirectMethod',true,'FinalXOR','0x00');
        case 10
            gen = comm.CRCGenerator('Polynomial',[1 1 0 0 1 0 0 1 1 1 1],'InitialConditions',zeros(1, 10),'DirectMethod',true,'FinalXOR',zeros(1, 10));
            det = comm.CRCDetector('Polynomial',[1 1 0 0 1 0 0 1 1 1 1],'InitialConditions',zeros(1, 10),'DirectMethod',true,'FinalXOR',zeros(1, 10));
        case 12
            gen = comm.CRCGenerator('Polynomial',[1 1 0 0 0 0 0 0 0 1 1 0 1],'InitialConditions',zeros(1, 12),'DirectMethod',true,'FinalXOR',zeros(1, 12));
            det = comm.CRCDetector('Polynomial',[1 1 0 0 0 0 0 0 0 1 1 0 1],'InitialConditions',zeros(1, 12),'DirectMethod',true,'FinalXOR',zeros(1, 12));
        case 16
            gen = comm.CRCGenerator('Polynomial',[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1],'InitialConditions',[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],'DirectMethod',true,...
                'FinalXOR',[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
            det = comm.CRCDetector('Polynomial',[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1],'InitialConditions',[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],'DirectMethod',true,...
                'FinalXOR',[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
        case 24
            gen = comm.CRCGenerator('Polynomial',[1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1],...
                'InitialConditions',[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],'DirectMethod',true,...
                'FinalXOR',[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
            det = comm.CRCDetector('Polynomial',[1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1],...
                'InitialConditions',[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],'DirectMethod',true,...
                'FinalXOR',[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
        otherwise
            error('Unsupported CRC length. Program terminates');
    end
end
