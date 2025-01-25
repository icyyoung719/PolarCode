function [gen, det] = get_crc_objective(crc_length)
% 根据指定的 CRC 长度创建 CRC 生成器和检测器对象
%   crc_length - CRC 校验的长度（支持 4, 6, 8, 10, 12, 16, 24）
%   gen - CRC 生成器对象，用于为数据添加 CRC 校验码
%   det - CRC 检测器对象，用于检测数据的 CRC 校验码
    switch crc_length
        case 4
            Polynomial = "x^4 + x + 1";
        case 6
            Polynomial = "x^6 + x + 1";
        case 8
            Polynomial = "x^8 + x^2 + x + 1";
        case 10
            Polynomial = "x^10 + x^9 + x^5 + x^4 + x + 1";
        case 12
            Polynomial = "x^12 + x^11 + x^3 + x^2 + 1";
        case 16
            Polynomial = "x^12 + x^11 + x^3 + x^2 + 1";
        case 24
            Polynomial = "x^24 + x^23 + x^6 + x^5 + x + 1";
        otherwise
            disp('Unsupported CRC length. Program terminates')
    end

     gen = comm.CRCGenerator(Polynomial = Polynomial);
     det = comm.CRCDetector(Polynomial = Polynomial);
end
