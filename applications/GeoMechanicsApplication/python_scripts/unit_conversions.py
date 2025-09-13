number_of_seconds_per_day = 86400.0


def seconds_to_days(duration_in_seconds):
    return duration_in_seconds / number_of_seconds_per_day


def days_to_seconds(duration_in_days):
    return duration_in_days * number_of_seconds_per_day


def Pa_to_kPa(stress_in_Pa):
    return stress_in_Pa / 1000.0
