def seconds_to_days(duration_in_seconds):
    return duration_in_seconds / (24.0 * 60.0 * 60.0)


def days_to_seconds(duration_in_days):
    return duration_in_days * 24.0 * 60.0 * 60.0


def Pa_to_kPa(stress_in_Pa):
    return stress_in_Pa / 1000.0