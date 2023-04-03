import os
import numpy as np
import pypulseq as pp
import yaml
from matplotlib import pyplot as plt

from mr_seq.pulse_seq import PulseSeq
from mr_seq.sequence import Sequence

from mr_seq.physical_setup import PhysicalSetup


class Sequence_gre:
    def check(self):
        return "Sequence_gre: checked!"

    def generate(
        self,
        matrix_shape: tuple,
        title: str,
        outdir: str,
        physics_configuration_filename: str,
        timestep: int = -1,
        verbal: bool = False,
    ):
        if verbal:
            print(
                f"Generate GRE sequence:\n\tnamed: '{title}' in directory: '{outdir}'"
            )

        #

        # Naming the GRE sequence output files for pypulseq and mri-sim format
        sequence_name = f"{title}_gre"
        sequence_pypulseq_filename = f"{sequence_name}.seq"
        sequence_config_filename = f"{sequence_name}.yml"

        # Create waveforms directory if not exists
        waveform_dir = os.path.join(outdir, "waveforms")
        os.makedirs(waveform_dir, exist_ok=True)

        # Physics system
        ps = PhysicalSetup(physics_configuration_filename)

        # Set raster time to timestep to ensure everything is with the same clock
        # TODO handle different raster times for rf and grad
        # raster_time: float = min(ps.grad_raster_time, ps.rf_raster_time)
        raster_time = timestep
        if timestep < 0:
            raster_time = ps.timestep

        # Init the system with the physics configuration
        system = pp.Opts(
            adc_dead_time=ps.adc_dead_time,
            gamma=ps.gamma,
            grad_raster_time=raster_time,
            grad_unit=ps.grad_unit,
            max_grad=ps.max_grad,
            max_slew=ps.max_slew,
            rf_dead_time=ps.rf_dead_time,
            rf_raster_time=raster_time,
            rf_ringdown_time=ps.rf_ringdown_time,
            rise_time=ps.rise_time,
            slew_unit=ps.slew_unit,
        )

        # Matrix size
        nz, ny, nx = matrix_shape

        # Define FOV and resolution
        fov = 200e-3
        alpha = 10  # alpha=°, flip_angle=rad
        slice_thickness = 3e-3

        # Timing
        TE = 4.3e-3  # np.array([4.3e-3])
        TR = 10e-3

        # RF spoiling increment (117 or 123)
        rf_spoiling_inc = 117

        # Create sequence in Pypulseq and mr-sim format
        seq_pypulseq = pp.Sequence(system=system)
        seq_yaml = Sequence(sequence_name)

        # TODO check TR and nTR exact meaning in simulator format
        seq_yaml.repetition_time = TR * ny
        seq_yaml.number_repetitions = 1

        # RF sinc pulse with slice select and slice select rephasing
        rf, gz, gzr = pp.make_sinc_pulse(
            flip_angle=alpha * np.pi / 180,
            duration=3e-3,  # TODO check ms
            slice_thickness=slice_thickness,
            apodization=0.5,
            time_bw_product=4,
            system=system,
            return_gz=True,
        )  # TBW=2 (rapid imaging), =4 (180), =8 (90), =12 (slab and saturation)
        # https://inst.eecs.berkeley.edu/~ee225e/sp13/notes/Lecture13.pdf

        # TODO pypulseq rf use 1e-6 hard coded for vector length
        rf_hard_code_clock = 1e-6
        rf_correction_ratio = int(raster_time / rf_hard_code_clock)
        rf_corrected_signal = rf.signal[0:-1:rf_correction_ratio].copy()
        rf_corrected_time = rf.t[: len(rf_corrected_signal)].copy()

        dur_rf = pp.calc_duration(rf)
        dur_gz = pp.calc_duration(gz)

        rf.signal = rf_corrected_signal
        rf.t = rf_corrected_time

        dur_rf = pp.calc_duration(rf)
        dur_gz = pp.calc_duration(gz)

        # Convert RF and associated gradient to waveforms of the same size in the format
        # of the GammaMRI-Simulator
        rf_waveform, rf_gz_waveform = Waveform.get_pulse_and_gradient_waveforms(
            rf,
            gz,
            os.path.join(waveform_dir, "rf.yml"),
            raster_time,
            self.PRECISION,
            ps.gamma,
        )

        # Define other gradients and ADC events
        delta_k = 1 / fov

        # Readout gradient
        gx = pp.make_trapezoid(
            channel="x", flat_area=nx * delta_k, flat_time=3.2e-3, system=system
        )

        # Convert gradient to waveform in the format of the GammaMRI-Simulator
        gx_waveform = TrapGradientWaveform(
            gx,
            os.path.join(waveform_dir, "gx.yml"),
            raster_time,
            self.PRECISION,
            ps.gamma,
        )

        adc = pp.make_adc(
            num_samples=nx, duration=gx.flat_time, delay=gx.rise_time, system=system
        )  # TODO useful to include in simulation ?

        # pre gradient in x
        gx_pre = pp.make_trapezoid(
            channel="x", area=-gx.area / 2, duration=1e-3, system=system
        )

        # Convert gradient to waveform in the format of the GammaMRI-Simulator
        gx_pre_waveform = TrapGradientWaveform(
            gx_pre,
            os.path.join(waveform_dir, "gx_pre.yml"),
            raster_time,
            self.PRECISION,
            ps.gamma,
        )

        # Rephasing gradient
        gz_reph = pp.make_trapezoid(
            channel="z", area=-gz.area / 2, duration=1e-3, system=system
        )

        # Convert gradient to waveform in the format of the GammaMRI-Simulator
        gz_reph_waveform = TrapGradientWaveform(
            gz_reph,
            os.path.join(waveform_dir, "gz_reph.yml"),
            raster_time,
            self.PRECISION,
            ps.gamma,
        )
        phase_areas = (np.arange(ny) - ny / 2) * delta_k

        # Define gradient spoiling
        gx_spoil = pp.make_trapezoid(channel="x", area=2 * nx * delta_k, system=system)

        # Convert gradient to waveform in the format of the GammaMRI-Simulator
        gx_spoil_waveform = TrapGradientWaveform(
            gx_spoil,
            os.path.join(waveform_dir, "gx_spoil.yml"),
            raster_time,
            self.PRECISION,
            ps.gamma,
        )

        gz_spoil = pp.make_trapezoid(
            channel="z", area=4 / slice_thickness, system=system
        )

        # Convert gradient to waveform in the format of the GammaMRI-Simulator
        gz_spoil_waveform = TrapGradientWaveform(
            gz_spoil,
            os.path.join(waveform_dir, "gz_spoil.yml"),
            raster_time,
            self.PRECISION,
            ps.gamma,
        )

        # Calculate timing
        delay_TE = (
            np.ceil(
                (
                    TE
                    - pp.calc_duration(gx_pre)
                    - gz.fall_time
                    - gz.flat_time / 2
                    - pp.calc_duration(gx) / 2
                )
                / seq_pypulseq.grad_raster_time
            )
            * seq_pypulseq.grad_raster_time
        )
        delay_TR = (
            np.ceil(
                (
                    TR
                    - pp.calc_duration(gz)
                    - pp.calc_duration(gx_pre)
                    - pp.calc_duration(gx)
                    - delay_TE
                )
                / seq_pypulseq.grad_raster_time
            )
            * seq_pypulseq.grad_raster_time
        )

        assert np.all(delay_TE >= 0)
        assert np.all(delay_TR >= pp.calc_duration(gx_spoil, gz_spoil))

        rf_phase = 0
        rf_inc = 0

        # CONSTRUCT SEQUENCE

        # Loop over phase encodes and define sequence blocks
        j = 0
        current_time: float = 0
        for i in range(ny):
            rf.phase_offset = rf_phase / 180 * np.pi
            adc.phase_offset = rf_phase / 180 * np.pi
            rf_inc = divmod(rf_inc + rf_spoiling_inc, 360.0)[1]
            rf_phase = divmod(rf_phase + rf_inc, 360.0)[1]

            # Pypulseq block
            seq_pypulseq.add_block(rf, gz)

            # Yaml equivalent
            pulse_seq = PulseSeq(current_time, rf_waveform.duration())
            pulse_seq.set_flip_angle(alpha)
            pulse_seq.set_b1(rf_waveform.filename)
            pulse_seq.set_gradient_z(rf_gz_waveform)
            seq_yaml.add_pulse_seq(pulse_seq)

            # Update time
            current_time += rf_waveform.duration()

            # Phase encoding gradient
            gy_pre = pp.make_trapezoid(
                channel="y",
                area=phase_areas[i],
                duration=pp.calc_duration(gx_pre),
                system=system,
            )

            # Convert gradient to waveform in the format of the GammaMRI-Simulator
            gy_pre_waveform = TrapGradientWaveform(
                gy_pre,
                os.path.join(waveform_dir, f"gy_pre{i}.yml"),
                raster_time,
                self.PRECISION,
                ps.gamma,
            )

            # Pypulseq block
            seq_pypulseq.add_block(gx_pre, gy_pre, gz_reph)

            # Yaml equivalent
            pulse_seq = PulseSeq(current_time, gx_pre_waveform.duration())
            pulse_seq.set_gradient_x(gx_pre_waveform)
            pulse_seq.set_gradient_y(gy_pre_waveform)
            pulse_seq.set_gradient_z(gz_reph_waveform)
            seq_yaml.add_pulse_seq(pulse_seq)

            # Update time
            current_time += gx_pre_waveform.duration()

            # Add delay TE
            seq_pypulseq.add_block(pp.make_delay(delay_TE))  # delay_TE[j]))
            current_time += float(delay_TE)

            # Pypulseq block
            seq_pypulseq.add_block(gx, adc)  # ADC is ignored for simulation

            # Yaml equivalent
            pulse_seq = PulseSeq(current_time, gx_waveform.duration())
            pulse_seq.set_gradient_x(gx_waveform)
            seq_yaml.add_pulse_seq(pulse_seq)

            # Update time
            current_time += gx_waveform.duration()

            gy_pre.amplitude = -gy_pre.amplitude

            # Convert gradient to waveform in the format of the GammaMRI-Simulator
            gy_spoil_waveform = TrapGradientWaveform(
                gy_pre,
                os.path.join(waveform_dir, f"gy_reph{i}.yml"),
                raster_time,
                self.PRECISION,
                ps.gamma,
            )

            # Pypulseq block
            seq_pypulseq.add_block(
                pp.make_delay(delay_TR), gx_spoil, gy_pre, gz_spoil
            )  # delay_TR[j]

            # Yaml equivalent : 3 events because 3 different durations
            pulse_seq = PulseSeq(current_time, gx_spoil_waveform.duration())
            pulse_seq.set_gradient_x(gx_spoil_waveform)
            seq_yaml.add_pulse_seq(pulse_seq)

            pulse_seq = PulseSeq(current_time, gy_spoil_waveform.duration())
            pulse_seq.set_gradient_y(gy_spoil_waveform)
            seq_yaml.add_pulse_seq(pulse_seq)

            pulse_seq = PulseSeq(current_time, gz_spoil_waveform.duration())
            pulse_seq.set_gradient_z(gz_spoil_waveform)
            seq_yaml.add_pulse_seq(pulse_seq)

            # Add delay to match TR
            current_time += float(delay_TR)

        (
            ok,
            error_report,
        ) = (
            seq_pypulseq.check_timing()
        )  # Check whether the timing of the sequence is correct
        if ok:
            print("Timing check passed successfully")
        else:
            print("Timing check failed. Error listing follows:")
            [print(e) for e in error_report]

        # ======
        # VISUALIZATION
        # ======
        seq_pypulseq.plot()

        # Trajectory calculation and plotting
        (
            ktraj_adc,
            ktraj,
            t_excitation,
            t_refocusing,
            t_adc,
        ) = seq_pypulseq.calculate_kspace()
        time_axis = np.arange(1, ktraj.shape[1] + 1) * system.grad_raster_time
        plt.figure()
        plt.plot(time_axis, ktraj.T)  # Plot the entire k-space trajectory
        plt.plot(t_adc, ktraj_adc[0], ".")  # Plot sampling points on the kx-axis
        plt.figure()
        plt.plot(ktraj[0], ktraj[1], "b")  # 2D plot
        plt.axis("equal")  # Enforce aspect ratio for the correct trajectory display
        plt.plot(ktraj_adc[0], ktraj_adc[1], "r.")  # Plot  sampling points
        plt.show()

        # Prepare the sequence output for the scanner
        seq_pypulseq.set_definition("FOV", [fov, fov, slice_thickness])
        seq_pypulseq.set_definition("Name", "gre")

        # Save .seq
        seq_pypulseq.write(os.path.join(outdir, sequence_pypulseq_filename))

        # Save yaml
        with open(os.path.join(outdir, sequence_config_filename), "w") as yaml_file:
            yaml.dump(seq_yaml.dict(), yaml_file, sort_keys=False)

        # Very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within
        # slew-rate limits
        # rep = seq_pypulseq.test_report()
        # print(rep)

        return sequence_config_filename
        #
